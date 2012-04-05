{-# LANGUAGE BangPatterns,ScopedTypeVariables #-}
import System.Environment (getArgs)
import System.Console.GetOpt
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
            ]

data Options = Options  {
        optAln  :: String,
        optTree :: String,
        optNumCats :: Int,
        optBootCount :: Int,
        optSeed :: Maybe Int,
        optLevel :: OptLevel
} deriving Show

defaultOptions :: Options
defaultOptions = Options {
        optAln = "in.fa",
        optTree = "in.tre",
        optNumCats = 4,
        optBootCount = 10,
        optSeed = Nothing,
        optLevel = BranchOpt
}

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
                  optLevel = optBoot
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
                (Just a,Right t,params)->   do let alpha:sigma:priorZero:[] = map read $ take 3 params
                                               let pAln = pAlignment a
                                               let (nSite,nCols,multiplicites) = getAlnData pAln where
                                                   getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts)
                                               let piF = fromList $ safeScaledAAFrequencies a
                                               let model = thmmModel (cats+1) wagS piF [priorZero,alpha,sigma]
                                               let t2 = addModelFx (structDataN (cats+1) AminoAcid (pAln) t) model [1.0]
                                               putStrLn "START"
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
                                               let pBEStr = getPartialBranchEnds t2
                                               let partials = map (map (map (toLists . trans))) $ toPBEList pBEStr
                                               let branchLengths = toPBELengths pBEStr
                                               let nBranch =  length $ toPBELengths pBEStr
                                               let sitelikes = [likelihoods t2] {-- FIXME, outer dimension is nProc, assume 1 --}
                                               let sitemap = mapBack pAln
                                               putStrLn $ show $ zip [0..] $ getCols pAln 
                                               let stochmapF a b c d e f g lM = calculateAndWrite nSite nState nBranch nProc nCols lM multiplicites a b c d e f g Nothing
                                               let stochmapT lM tree = stochmapF sitemap' partials' qset' sitelikes' pi_i' branchLengths' mixProbs' lM where
                                                                    sitemap' = mapBack $ pAlignment pAln'
                                                                    pAln' = getAln tree
                                                                    pBEStr' = getPartialBranchEnds tree
                                                                    partials' = map (map (map (toLists . trans))) $ toPBEList pBEStr'
                                                                    qset' = map toLists $ getQMat tree
                                                                    sitelikes' = [likelihoods tree] --FIXME
                                                                    pi_i' = map toList $ getPi tree
                                                                    branchLengths' = toPBELengths pBEStr
                                                                    mixProbs' = getPriors tree
                                               let numQuantile = 500
                                               let stdGens = take numSim $ genList stdGen
                                               let alnLength = length $ Phylo.Likelihood.columns pAln
                                               let simulations' = map (\x->makeSimulatedTree AminoAcid (cats+1) x alnLength t2) stdGens
                                               simulations :: [DNode] <- case optBoot of 
                                                                      FullOpt -> do mVar <- newEmptyMVar
                                                                                    putStrLn "FULL OPT"
                                                                                    let threads = Sync.numCapabilities
                                                                                    stagger <- replicateM threads newEmptyMVar
                                                                                    mapM_ (\x->putMVar x ()) stagger
                                                                                    let staggerVars = concat $ repeat stagger
                                                                                    let opt mv (mytree,hv,id) = do gottheball<-takeMVar hv
                                                                                                                   startTime <- getCurrentTime
                                                                                                                   putStrLn $ "Started " ++ (show id) ++ " " ++ (show startTime)
                                                                                                                   hFlush stdout
                                                                                                                   ans <- optBSParamsBLIO (1,numModels) (makeMapping (allIn []) mytree) (map (\x->0.01) lower) lower upper [1.0] myModel mytree ([sigma,priorZero,alpha])
                                                                                                                   putMVar hv ()
                                                                                                                   putMVar mv $ fst ans 
                                                                                                                   endTime <- getCurrentTime
                                                                                                                   let diffTime = endTime `diffUTCTime` startTime
                                                                                                                   putStrLn $ "Time taken " ++ (show (diffTime)) 
                                                                                                                   hFlush stdout where
                                                                                                                     lower = (replicate numModels $ Just 0.0) ++ [Just 0.001,Just 0.001]                                                                                                           
                                                                                                                     upper = (replicate numModels Nothing) ++ [Just 0.99,Nothing]                                                                                                                  
                                                                                                                     numModels = 1
                                                                                                                     myModel = thmmPerBranchModel (cats+1) wagS piF
                                                                                    mapM_ (forkIO . opt mVar) $ zip3 simulations' staggerVars [0..]
                                                                                    replicateM (length simulations') $ takeMVar mVar 
                                                                      BranchOpt -> do putStrLn "BRANCH OPT"
                                                                                      return $ parMap rseq optBLDFull0 simulations'
                                                                      QuickBranchOpt -> do putStrLn "QUICK BRANCH OPT (todo)"
                                                                                           return $ parMap rseq optBLDFull0 simulations'
                                                                      NoOpt -> return $ simulations'
                                               let outputMat (name,lM) = do ans' <- stochmapT lM t2 -- sitemap partials qset sitelikes pi_i branchLengths mixProbs Nothing
                                                                            let ans = stochmapOrder ans' sitemap (getPriors t2) 
                                                                            simStoch' <- mapM (stochmapT lM) simulations
                                                                            let simStoch = map (\(x,y)->stochmapOrder x (treeToMap y) (getPriors y)) $ zip simStoch' simulations
                                                                            --print $ map fst simStoch
                                                                            let simStochDesc = map discrep simStoch
                                                                            let ansDesc = discrep ans
                                                                            let tot x = foldr (+) (0.0) x
                                                                            let (line,lower,upper) = makeQQLine numQuantile (map concat simStochDesc) (concat ansDesc)
                                                                            let (line1,lower1,upper1) = makeQQLine numQuantile (map (map tot) simStochDesc) (map tot ansDesc)
                                                                            let (line2,lower2,upper2) = makeQQLine numQuantile (map (map tot)  $ map transpose simStochDesc) (map tot $ transpose ansDesc)
                                                                            let edgeQuantileMap = zip (map (\(a,b,c,d,e) -> d) $ getPartialBranchEnds t2) (perLocationQuantile numQuantile (map (map tot) simStochDesc) (map tot ansDesc))
                                                                            --print edgeQuantileMap
                                                                            let tempTree = annotateTreeWith edgeQuantileMap t2
                                                                            writeFile (name ++ "-colour-tree.xml")  $ unlines $ quantilePhyloXML (annotateTreeWith edgeQuantileMap t2)
                                                                            renderableToPDFFile (makePlot line (zip lower upper) PDF) 480 480 $ name ++ "-out-all.pdf"
                                                                            renderableToPDFFile (makePlot line1 (zip lower1 upper1) PDF) 480 480 $ name ++ "-out-branch.pdf"
                                                                            renderableToPDFFile (makePlot line2 (zip lower2 upper2) PDF) 480 480 $ name ++ "-out-site.pdf"
                                               outputMat ("inter",interLMat (cats+1) 20)
                                               outputMat ("intra",intraLMat (cats+1) 20)
                                               return Nothing
          case output of
               Just str -> putStrLn str
               Nothing -> return ()

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


                  
makeQQLine :: Int -> [[Double]] -> [Double] -> ([Double],[Double],[Double])
makeQQLine numQuantile simulated empirical = (revQuantile empirical,lower,upper) where
                                        revQuantile x = map (quantileF (revQuantileRawF x)) [0..numQuantile]
                                        revQuantileRawF x = UVec.fromList $ map (\y->(fromIntegral $ reverseQuantile simQuantiles y)/(fromIntegral numQuantile)) x

                                        quantileF d q = continuousBy medianUnbiased q numQuantile d
                                        simVec = UVec.fromList $ concat simulated
                                        simQuantiles = map (quantileF simVec) [0..numQuantile]

                                        numSim = length simulated
                                        revQuantileSimS = map sort $ transpose $ map revQuantile simulated
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
optParamsAndBL model tree params priors lower upper cutoff = optWithBS' [] cutoff (0,0) Nothing lower upper priors model (dummyTree tree) params                                                            
optParamsAndBLIO model tree params priors lower upper cutoff = optWithBSIO' [] cutoff (0,0) Nothing (map (\x->0.01) lower) (map (\x->1E-4) lower) lower upper priors model (dummyTree tree) params                                                            

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
                                                                                                                                                                                                                                              
                                                             
