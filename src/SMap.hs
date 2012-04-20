{-# LANGUAGE BangPatterns,ScopedTypeVariables #-}
import System.Environment (getArgs)
import System.Console.GetOpt
import System.Exit
import Phylo.Alignment
import Data.List
import Control.Monad
import Control.Concurrent
import qualified GHC.Conc.Sync as Sync
import Phylo.Tree
import Phylo.Data
import Phylo.Likelihood
import Phylo.Opt
import Numeric.GSL.Minimization
import Data.Packed.Vector
import Debug.Trace
import System.Random
import Data.Packed.Matrix
import System.IO
import System.Posix.Types
import System.Posix.IO
import Foreign.C.String (CString,newCString)
import Phylo.Matrix
import Stochmap
import Data.Char (isSpace)
import Phylo.PhyloXML
import Statistics.Quantile
import Phylo.Graphics.Plotting
import Graphics.Rendering.Chart.Renderable (renderableToPDFFile)
import Control.DeepSeq
import Data.Char
import Control.Parallel.Strategies
import Data.Time.Clock
import Phylo.NLOpt
import qualified Data.Vector.Unboxed as UVec

data LMat = Inter | Intra deriving Show
data OptLevel = FullOpt | BranchOpt | QuickBranchOpt | NoOpt deriving Show

--handle options like http://leiffrenzel.de/papers/commandline-options-in-haskell.html
options = [ Option ['a'] ["alignment"] (ReqArg optAlnR "FILE") "Alignment"
          , Option ['t'] ["tree"] (ReqArg optTreeR "FILE") "Tree"
          , Option [] ["num-cats"] (ReqArg optNumCatsR "NUMGAMMACATS") "Number of Gamma Rate Categories" 
          , Option ['n'] ["bootstrap"] (ReqArg optBootCountR "BOOTSTRAPS") "Number of bootstraps to perform"
          , Option ['s'] ["seed"] (ReqArg optSeedR "SEED") "RNG seed"
          , Option [] ["opt-bootstrap"] (ReqArg optFR "full, branch, quick, none") "Optimisation of bootstraps"
          , Option [] ["thmm"] (OptArg thmmModelR "init params") "Use THMM model"
          , Option [] ["opt"] (OptArg optAlgR "opt function") "Optimise model for real data"
            ]

thmmModelR arg opt | trace "THMM found" True  = case arg of 
                      Nothing | trace "THMM NOTHING" True -> return opt {optModel = Thmm 1.0 0.1 1.0 }
                      Just args | trace ("THMM " ++ args) True -> case (map read $ splitBy ' ' args) of 
                                        [a,b,c] -> return opt { optModel = Thmm a b c }
                                        _       -> error $ "Can't parse three doubles from " ++ args

optAlgR arg opt = case arg of 
                        Nothing -> return opt {optAlg = OptMethod var2}
                        Just name -> return opt {optAlg = OptMethod met} where 
                                        met = case (map toLower name) of 
                                              "bobyqa" -> bobyqa
                                              "cobyla" -> cobyla
                                              "mma" -> mma
                                              "slsqp" -> slsqp
                                              "newton" -> newton
                                              "var1" -> var1
                                              "var2" -> var2
                                              "lbfgs" -> lbfgs
                                              "sbplx" -> sbplx
                                              "neldermead" -> neldermead
                                              "praxis" -> praxis
                                              "newuoa" -> newuoa
                                              _ -> error $ "Can't parse optimisation method " ++ name

data Options = Options  {
        optAln  :: String,
        optTree :: String,
        optNumCats :: Int,
        optBootCount :: Int,
        optSeed :: Maybe Int,
        optLevel :: OptLevel,
        optModel :: Model,
        optAlg :: OptAlg
} deriving Show

defaultOptions :: Options
defaultOptions = Options {
        optAln = "in.fa",
        optTree = "in.tre",
        optNumCats = 4,
        optBootCount = 10,
        optSeed = Nothing,
        optLevel = BranchOpt,
        optModel = Thmm 1.0 0.1 1.0,
        optAlg = OptNone
}

-- alpha sigma pInv | alpha
data Model = Thmm Double Double Double | Ras Double deriving Show
data OptAlg = OptNone | OptMethod NLOptMethod deriving Show

optAlnR arg opt = return opt { optAln = arg }
optTreeR arg opt = return opt { optTree = arg }
optNumCatsR arg opt = return opt {optNumCats = (read arg)}
optBootCountR arg opt = return opt {optBootCount = (read arg)}
optSeedR arg opt = return opt {optSeed = Just (read arg)}
optFR arg opt = optFR' (map toLower arg) opt
optFR' arg opt | arg=="full" = return opt {optLevel = FullOpt}
               | arg=="branch" = return opt {optLevel = BranchOpt}
               | arg=="quick" = return opt {optLevel = QuickBranchOpt}
               | arg=="none" = return opt {optLevel = NoOpt}
               | otherwise = error $ "need to specify one of full, branch, quick or none instead of " ++ arg

---

main = do args <- getArgs
          let ( actions, nonOpts, msgs ) = getOpt Permute options args
          opts <- foldl (>>=) (return defaultOptions) actions
          print opts
          let Options {
                  optAln = aln,
                  optTree = tree,
                  optNumCats = cats,
                  optBootCount = numSim,
                  optSeed = seed,
                  optLevel = optBoot,
                  optModel = modelParams,
                  optAlg = optMethod
          } = opts
          print seed
          print numSim
          print cats 
          stdGen <- case seed of
                    Nothing -> getStdGen
                    Just x -> return $ mkStdGen x
          aln' <- parseAlignmentFile parseUniversal aln                                                                                                       
          tree' <- (liftM readBiNewickTree) (readFile tree)                                                                                                    
          output <- case (aln',tree',nonOpts) of 
                (Nothing,_,_) -> do putStrLn "Failed to parse alignment"
                                    return Nothing
                (_,Left err,_) -> do putStrLn $ "Failed to parse tree " ++ err
                                     return Nothing
                (Just a,Right t,params)->   do let (sigma',priorZero',alpha') = case modelParams of 
                                                                                 Thmm a b c -> (a,b,c)
                                               let pAln = pAlignment a
                                               let piF = fromList $ safeScaledAAFrequencies a
                                               let model = thmmModel (cats+1) wagS piF [priorZero',alpha',sigma']
                                               let t2' = addModelFx (structDataN (cats+1) AminoAcid (pAln) t) model [1.0]
                                               let (t2,[priorZero,alpha,sigma]) = case optMethod of 
                                                                                       OptNone                        -> (t2',[priorZero',alpha',sigma'])
                                                                                       (OptMethod method)             -> optThmmModel 1 cats piF wagS t2' priorZero' alpha' sigma' method 
                                               print t2
                                               let a = map (fst . leftSplit) $ getAllF t2
                                               let b = getLeftSplit t2
                                               print $ "OK? " ++ (show (a==b))
                                               putStrLn "END"
                                               let nState = (cats+1)*20
                                               let nProc = 1 {-- fixme --}
                                               let pi_i = map toList (snd model) 
                                               let mixProbs = getPriors t2
                                               let qMat = thmmModelQ (cats+1) wagS piF [priorZero,alpha,sigma]
                                               let qset = [toLists qMat]
                                               let sitelikes = [likelihoods t2] {-- FIXME, outer dimension is nProc, assume 1 --}
                                               let sitemap = mapBack pAln
                                               putStrLn $ show $ zip [0..] $ getCols pAln 
                                               let stochmapF tree a b c d e f g lM = calculateStochmap nSite nState nBranch nProc nCols lM multiplicites a b c d e f g where
                                                                                (nSite,nCols,multiplicites) = getAlnData pAln
                                                                                pAln = pAlignment $ getAln tree
                                                                                pBEStr = getPartialBranchEnds tree
                                                                                nBranch =  length $ toPBELengths pBEStr
                                               let stochmapT lM tree = stochmapF tree sitemap' partials' qset' sitelikes' pi_i' branchLengths' mixProbs' lM where
                                                                    sitemap' = mapBack $ pAlignment pAln'
                                                                    pAln' = getAln tree
                                                                    pBEStr' = getPartialBranchEnds tree
                                                                    partials' = map (map (map (toLists . trans))) $ toPBEList pBEStr'
                                                                    qset' = map toLists $ getQMat tree
                                                                    sitelikes' = [likelihoods tree] --FIXME
                                                                    pi_i' = map toList $ getPi tree
                                                                    branchLengths' = toPBELengths pBEStr'
                                                                    mixProbs' = getPriors tree
                                               let stochmapTT lM tree = discrep $ stochmapOrder (stochmapT lM tree) (treeToMap tree) (getPriors tree)
                                               let numQuantile = 500
                                               let stdGens = take numSim $ genList stdGen
                                               let alnLength = length $ Phylo.Likelihood.columns pAln
                                               let simulations' = map (\x-> patternSimulation t t2 model [1.0] x AminoAcid cats alnLength) stdGens
                                               let simulations :: [DNode] = case (optBoot,optMethod) of 
                                                                      (FullOpt,OptMethod method) -> map (\x->fst $ optThmmModel 1 cats piF wagS x priorZero alpha sigma method) simulations'
                                                                      (FullOpt,OptNone) -> map (\x->fst $ optThmmModel 1 cats piF wagS x priorZero alpha sigma bobyqa) simulations'
                                                                      (BranchOpt,_) -> map optBLDFull0 simulations'
                                                                      (QuickBranchOpt,_) -> map optBLDFull0 simulations'
                                                                      (NoOpt,_) -> simulations'
                                               let simS = map (dualStochMap stochmapTT (interLMat (cats+1) 20) (intraLMat (cats+1) 20))  simulations
                                               let (simStochInter,simStochIntra) = unzip (simS `using` (parListChunk (numSim `div` Sync.numCapabilities) rdeepseq))
                --                               let (simStochInter,simStochIntra) = unzip simS --using` (parListChunk (numSim `div` Sync.numCapabilities) rdeepseq))
                                               let ansInter = stochmapTT (interLMat (cats+1) 20) t2
                                               let ansIntra = stochmapTT (intraLMat (cats+1) 20) t2
                                               let outputMat (name,simStochDesc,ansDesc) = do let tot x = foldr (+) (0.0) x
                                                                                              let (line,lower,upper,pval) = makeQQLine numQuantile (map concat simStochDesc) (concat ansDesc)
                                                                                              let (line1,lower1,upper1,pval1) = makeQQLine numQuantile (map (map tot) simStochDesc) (map tot ansDesc)
                                                                                              let (line2,lower2,upper2,pval2) = makeQQLine numQuantile (map (map tot)  $ map transpose simStochDesc) (map tot $ transpose ansDesc)
                                                                                              let edgeQuantileMap = zip (map (\(a,b,c,d,e) -> d) $ getPartialBranchEnds t2) (perLocationQuantile numQuantile (map (map tot) simStochDesc) (map tot ansDesc))
                                                                                              let tempTree = annotateTreeWith edgeQuantileMap t2
                                                                                              writeFile (name ++ "-colour-tree.xml")  $ unlines $ quantilePhyloXML (annotateTreeWith edgeQuantileMap t2)
                                                                                              writeFile (name ++ "-pvals.txt")  $ unlines $ map (\(x,y) -> x ++ " " ++ (show y)) $ zip ["all","edge","site"] [pval,pval1,pval2]
                                                                                              renderableToPDFFile (makePlot line (zip lower upper) PDF) 480 480 $ name ++ "-out-all.pdf"
                                                                                              renderableToPDFFile (makePlot line1 (zip lower1 upper1) PDF) 480 480 $ name ++ "-out-edge.pdf"
                                                                                              renderableToPDFFile (makePlot line2 (zip lower2 upper2) PDF) 480 480 $ name ++ "-out-site.pdf"
                                               outputMat ("switch",simStochInter,ansInter)
                                               outputMat ("subs",simStochIntra,ansIntra)
                                               return Nothing
          case output of
               Just str -> putStrLn str
               Nothing -> return ()



dualStochMap stochmapTT a b x = (stochmapTT a x, stochmapTT b x)
length2 = map length
length3 = map length2

--FIXME not sure about this
reverseQuantile myQList point = case (filter (\(x,y) -> point <= y) $ zip [0..] myQList) of
                                        (x:xs) -> fst $x
                                        [] -> length $ myQList

discrep :: ([[Double]], [Double]) -> [[Double]]
discrep ((x:xs),(y:ys)) = (map (\i->i-y) x):(discrep (xs,ys))
discrep ([],[]) = []

perLocationQuantile :: Int -> [[Double]] -> [Double] -> [Double]
perLocationQuantile numQuantile simulated locdist = map ((/ (fromIntegral numQuantile)) . fromIntegral . (reverseQuantile simQuantiles)) locdist where
        quantileF d q = continuousBy medianUnbiased q numQuantile d
        simVec = UVec.fromList $ concat simulated
        simQuantiles = map (quantileF simVec) [0..numQuantile]


                  
makeQQLine :: Int -> [[Double]] -> [Double] -> ([Double],[Double],[Double],Double)
makeQQLine numQuantile simulated empirical = (empQuantile,lower,upper,pvalue) where
                                        revQuantile x = map (quantileF (revQuantileRawF x)) [0..numQuantile]
                                        revQuantileRawF x = UVec.fromList $ map (\y->(fromIntegral $ reverseQuantile simQuantiles y)/(fromIntegral numQuantile)) x

                                        quantileF d q = continuousBy medianUnbiased q numQuantile d
                                        simVec = UVec.fromList $ concat simulated
                                        simQuantiles = map (quantileF simVec) [0..numQuantile]

                                        numSim = length simulated
                                        revQuantileSimS = map sort $ transpose revQuantileSimS'
                                        revQuantileSimS' = map revQuantile simulated
                                        empQuantile = revQuantile empirical
                                        pvalue = cramerVonMises empQuantile revQuantileSimS'

                                        lowerI = floor $ (fromIntegral numSim) * 0.025
                                        upperI = (ceiling ((fromIntegral numSim) * 0.975)) - 1 
                                        lower = map (!!lowerI) revQuantileSimS
                                        upper = map (!!upperI) revQuantileSimS

treeToMap tree = mapBack $ pAlignment $ getAln tree


splitBy delim s = ans where
  (token,rest) = span (/=delim) s
  ans = case rest of
    [] -> [token]
    t  -> token : splitBy delim (tail rest)

splitsStr = unlines . splitsStr' 0
splitsStr' i [] = []
splitsStr' i ((l,r):xs) = ((show i) ++ " " ++ (intercalate " " l) ++ " | " ++ (intercalate " " r)) : (splitsStr' (i+1) xs)

getCols (PatternAlignment _ _ c _ _) = c

normalise list = map ( / total) list where     
                 total = foldr (+) 0.0 list

safeScaledAAFrequencies = normalise . map (\x-> if x < 1e-15 then 1e-15 else x) . scaledAAFrequencies

trim = f . f where 
   f = reverse . dropWhile isSpace
clean = clean' . trim where
   clean' (' ':' ':xs) = clean' (' ':xs)
   clean' (x:xs) = x:(clean' xs)
   clean' [] = []
joinWith i (x:x':xs) = (show x) ++ i ++ (joinWith i (x':xs))
joinWith i [x] = (show x) 
joinWith i [] = []

stochmapOrder :: ([[[Double]]],[[Double]]) -> [Int] -> [Double] -> ([[Double]],[Double])
stochmapOrder (condE,priorE) mapping priors = order where
                                                condE' = map (fixProc2 priors) condE
                                                priorE' = map (fixProc priors) priorE
                                                fixProc pr x = (foldr (+) (0.0) (map (\(x,y) -> x*y) $ zip pr x))
                                                fixProc2 pr xs = (map $ fixProc pr) (transpose xs)
                                                order = (transpose $ map (\i-> (map (!!i) condE')) mapping, priorE')

stochmapOut :: ([[[Double]]],[[Double]]) -> [Int] -> [Double] -> (String -> IO()) -> IO ()
stochmapOut (condE',priorE') mapping priors f = do 
                                                f header
                                                (mapM . mapM) f remainder
                                                f "\n"
                                                f headerSite 
                                                mapM f remainderSite
                                                f "\n"
                                                f headerBranch
                                                mapM f remainderBranch
                                                return () where
                                                 (condE''',priorE''') = stochmapOrder (condE',priorE') mapping priors
                                                 priorE'' = replicate (length mapping) priorE'''
                                                 condE'' = transpose condE'''
                                                 condE=transpose condE''
                                                 priorE=transpose priorE''
                                                 header = "Branch\tSite\tConditional_expectation\tPrior_expectation\n"
                                                 remainder = fmtBranch [0..] $ zip condE priorE 
                                                 fmtBranch (b:bs) ((cond,prior):xs) = (fmtBranch' b [0..] cond prior)  : (fmtBranch bs xs)
                                                 fmtBranch bs [] = []
                                                 fmtBranch' b sites [] priors = []
                                                 fmtBranch' b (site:sites) (cond:cs) (prior:ps) = ((show b) ++ "\t" ++ (show site) ++ "\t" ++ (show cond) ++ "\t" ++ (show prior) ++ "\n") : (fmtBranch' b sites cs ps)
                                                 headerBranch = "Branch\tTotal_conditional_expectation\tTotal_prior_expectation\n"
                                                 remainderBranch = fmtBranch2 [0..] $ zip condE priorE
                                                 total xs = foldr (+) 0.0 xs
                                                 fmtBranch2 bs [] = []
                                                 fmtBranch2 (b:bs) ((cond,prior):xs) = ((show b) ++"\t" ++ (show $ total cond) ++ "\t" ++ (show $ total prior) ++ "\n") : (fmtBranch2 bs xs)
                                                 headerSite = "Site\tTotal_conditional_expectation\tTotal_prior_expectation\n"
                                                 remainderSite = fmtBranch2 [0..] $ zip condE'' priorE''

qList qFunc numQ = map (\x-> qFunc x numQ) [0..numQ]




uniformQ :: [([[Double]],[Double])] -> (Int -> Int -> Double)
uniformQ simData = partialFunc where 
                      partialFunc q k | trace ("Made vec " ++ (show vec)) True = continuousBy medianUnbiased q k vec
                      vec | trace ("Making vector " ++ (take 10 $ show simData)) True = (UVec.fromList $ concatMap allDisc simData)

linearFromRaw mapping priors ans | trace "linearFromRaw" True = linearAns $ stochmapOrder ans mapping priors
linearAns reformattedAns = allDisc reformattedAns
reformatAns mapping priors stochResult = (map $ uncurry3 stochmapOrder) $ zip3 stochResult mapping priors
uniformQRaw mapping priors stochResult = uniformQ $ reformatAns mapping priors stochResult

uncurry3 f (a,b,c) = f a b c 

allDisc :: ([[Double]],[Double]) -> [Double]
allDisc a | trace ("allDisc" ++ (take 20 $ show a)) False = undefined
allDisc (condE:xs,priorE:ys) | traceShow ("len " ++ (show $ 1 + (length xs)) ++ " " ++ (show $ 1 + (length ys))) True = (allDisc' condE priorE) ++ (allDisc (xs,ys))
allDisc ([],[]) = []
allDisc (c,p) | traceShow (c,p) True = undefined
--allDisc' condE priorE | trace ("prior cond " ++ (show priorE) ++ "  " ++ (show condE)) False = undefined
allDisc' condE priorE  = map (\x -> priorE - x) condE
                                        
allInGeneric method sets (l,r) = case (findIndex (==True) $ map (method l r) sets) of                                                                                                                                                         
                                        Just a -> a                                                                                                                                                                                           
                                        Nothing -> other                                                                                                                                                                                      
                                 where other = length sets                                                                                                                                                                                    

allIn = allInGeneric allInDynamic'                                                                                                                                                                                                            
allInNoRoot = allInGeneric allInNoRoot'                                                                                                                                                                                                       
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
allIn' :: [String] -> [String] -> [String] -> Bool                                                                                                                                                                                            
allIn' l r set = (l \\ set == [] || r \\ set == []) --we want l or r to be fully contained within 'set'                                                                                                                                       
                                                                                                                                                                                                                                              
allInNoRoot' :: [String] -> [String] -> [String] -> Bool                                                                                                                                                                                      
allInNoRoot' l r set = ((l \\ set == []) && (set \\ l /= [])) || ((r \\ set == []) &&  (set \\ r /= [])) --we want l or r to be fully contained with in 'set', but not the same as 'set'                                                      
                                                                                                                                                                                                                                              
--if first symbol is '@' then use noRoot version                                                                                                                                                                                              
allInDynamic' :: [String] -> [String] -> [String] -> Bool                                                                                                                                                                                     
allInDynamic' l r ("@":set) = allInNoRoot' l r set                                                                                                                                                                                            
allInDynamic' l r set = allIn' l r set                                                                                                                                                                                                        
                                                                                                                                                                                                                                              
                                                             
optThmmModel numModels cats pi s t2 priorZero alpha sigma method = optBSParamsBLIO method (0,numModels) (makeMapping (allIn []) t2) (map (\x->0.01) lower) lower upper [1.0] myModel t2 ([priorZero,alpha,sigma]) where
        lower = [Just 0.01,Just 0.1,Just 0.01]                                                                                                           
        upper = [Just 0.99,Just 200.0,Just 500.0]                                                                                                                  
        myModel = thmmModel (cats+1) s pi

--this is a bit of a hack really
--all the data needed is in modelTree
--but this hack will 'recompress' the data
--to site patterns
patternSimulation origTree modelTree model priors stdGen dataType cats len = trace ("cols " ++ (show len) ++ " " ++ (show $ length (Phylo.Alignment.columns lAln))) ans where
                             ans = addModelFx (structDataN (cats+1) AminoAcid (pAln) origTree) model priors 
                             pAln = pAlignment $ lAln
                             lAln = getAln $ makeSimulatedTree AminoAcid (cats+1) stdGen len modelTree 
getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts)
cramerVonMises empQ theQs = pvalue where
        v = fromIntegral $ length empQ
        cram x = foldr (+) 0.0 $ map (\(w,i)-> w - (fromIntegral ((2*i)-1)/(2*v)) + (1.0/(12*v))) $ zip x (map fromIntegral [1..])
        wv2_emp = cram empQ
        wv2_thes = map cram theQs
        pvalue = (fromIntegral $ length (filter (<=wv2_emp) wv2_thes)) / (fromIntegral $ length theQs)
