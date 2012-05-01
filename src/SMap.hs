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
import System.ProgressBar
import Text.Printf
import qualified Data.Vector.Unboxed as UVec

data LMat = Inter | Intra deriving Show
data OptLevel = FullOpt Double | BranchOpt | QuickBranchOpt | NoOpt deriving Show
data SubModel = WAG | WAGF | JTT | JTTF | CustomS String | CustomSF String deriving Show

--handle options like http://leiffrenzel.de/papers/commandline-options-in-haskell.html
options = [ Option ['a'] ["alignment"] (ReqArg optAlnR "FILE") "Alignment"
          , Option ['t'] ["tree"] (ReqArg optTreeR "FILE") "Tree"
          , Option [] ["num-cats"] (ReqArg optNumCatsR "NUMGAMMACATS") "Number of Gamma Rate Categories" 
          , Option ['n'] ["bootstrap"] (ReqArg optBootCountR "BOOTSTRAPS") "Number of bootstraps to perform"
          , Option ['s'] ["seed"] (ReqArg optSeedR "SEED") "RNG seed"
          , Option [] ["opt-bootstrap"] (ReqArg optFR "full, branch, quick, none") "Optimisation of bootstraps"
          , Option [] ["thmm"] (OptArg thmmModelR "init params") "Use THMM model"
          , Option [] ["gamma"] (OptArg gammaModelR "init params") "Use Gamma model"
          , Option [] ["opt"] (OptArg optAlgR "opt function") "Optimise model for real data (default)"
          , Option [] ["noopt"] (NoArg optNoneR) "Don't optimise model for real data"
          , Option [] ["raw"] (NoArg optRawR) "include raw data output"
          , Option ['D'] ["debug"] (NoArg optLogR) "debugging output"
          , Option ['h'] ["help"] (NoArg printHelp) "print help"
          , Option ['m'] ["model"] (ReqArg optSubR "SUBMODEL") "model (wag or wag,F (default) or jtt or jtt,F or <customFile> or <customFile>,F)"
            ]


printHelp opt = do putStrLn $ usageInfo header options
                   putStrLn example
--                   putStrLn reference
                   exitSuccess where
                       header = unlines ["SMap: Graphical checking of phylogenetic models using stochastic mapping"
                                        ,""]
                       example = unlines ["Examples:","Optimise and check a THMM model using 200 bootstrap replicates"
                                         ,"Using alignment ef1.phy and tree ef1.tre"
                                         ,""
                                         ,"smap -a ef1.phy -t ef1.tre --thmm --bootstrap 50 --opt-bootstrap --opt"
                                         ]
--                       reference = unlines ["Graphical checking of covarion models using stochastic mapping"
--                                           ,"Benjamin P. Blackburne, Simon Whelan, and Matthew Spencer"]
optRawR opt = return opt {optRaw = True}
optLogR opt = return opt {optLog = Logger $ hPutStrLn stderr,
                          optDebug = True}

optSubR arg opt = do let model = case arg of 
                                        "jtt" -> JTT
                                        "jtt,F" -> JTTF
                                        "wag" -> WAG
                                        "wag,F" -> WAGF
                                        x | (isF x) -> CustomSF (reverse $ drop 2 $ reverse $ x)
                                          | otherwise -> CustomS x
                     return opt {optSub = model}  where
                     isF (x:y:[]) = [x,y] == ",F"
                     isF (x:xs) = isF xs
                     isF [] = False

                                        
                                        
thmmModelR arg opt | trace "THMM found" True  = case arg of 
                      Nothing | trace "THMM NOTHING" True -> return opt {optModel = Thmm 0.1 1.0 1.0 }
                      Just args | trace ("THMM " ++ args) True -> case (map read $ splitBy ' ' args) of 
                                        [a,b,c] -> return opt { optModel = Thmm a b c }
                                        _       -> error $ "Can't parse three doubles from " ++ args

gammaModelR arg opt | trace "THMM found" True  = case arg of 
                      Nothing | trace "THMM NOTHING" True -> return opt {optModel = Ras 1.0 }
                      Just args | trace ("THMM " ++ args) True -> case (read args) of 
                                        a -> return opt { optModel = Ras a }


optNoneR opt = return opt {optAlg = OptNone}
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
        optAln  :: Maybe String,
        optTree :: Maybe String,
        optNumCats :: Int,
        optBootCount :: Int,
        optSeed :: Maybe Int,
        optLevel :: OptLevel,
        optModel :: Model,
        optAlg :: OptAlg,
        optRaw :: Bool,
        optLog :: Logger,
        optDebug :: Bool,
        optSub :: SubModel
} deriving Show

nullOut :: Logger
nullOut = Logger n where
        n _  = return ()

newtype Logger = Logger (String -> IO ())
instance Show Logger where
        show x = ""
defaultOptions :: Options
defaultOptions = Options {
        optAln = Nothing,
        optTree = Nothing,
        optNumCats = 4,
        optBootCount = 10,
        optSeed = Nothing,
        optLevel = FullOpt 1E-1,
        optModel = Thmm 0.1 1.0 1.0,
        optAlg = OptMethod var2,
        optRaw = False,
        optLog = nullOut,
        optDebug = False,
        optSub = WAGF
}

-- alpha sigma pInv | alpha
data Model = Thmm Double Double Double | Ras Double deriving Show
data OptAlg = OptNone | OptMethod NLOptMethod deriving Show

optAlnR arg opt = return opt { optAln = Just arg }
optTreeR arg opt = return opt { optTree = Just arg }
optNumCatsR arg opt = return opt {optNumCats = (read arg)}
optBootCountR arg opt = return opt {optBootCount = (read arg)}
optSeedR arg opt = return opt {optSeed = Just (read arg)}
optFR arg opt = optFR' (map toLower arg) opt
optFR' arg opt | (take 4 arg)=="full" = return opt {optLevel = FullOpt (read $ drop 4 arg)}
               | arg=="branch" = return opt {optLevel = BranchOpt}
               | arg=="quick" = return opt {optLevel = QuickBranchOpt}
               | arg=="none" = return opt {optLevel = NoOpt}
               | otherwise = error $ "need to specify one of full, branch, quick or none instead of " ++ arg

---

printNiceOpt len (x,x') = do
   let string = getNiceOpt2 len x x'
   let cr = string `deepseq` '\r'
   putChar cr
   putStr (replicate len ' ')
   putChar cr
   putStr string

getNiceOpt len (a,b,c) = outL ++ outR where
                                outL = (show c) ++ " : " ++ (format $ logLikelihood a)
                                outR = replicate (len - (length outL)) ' '
                                format = printf "%.5f"

getNiceOpt2 len (a,b,c) (a',b',c') = outL ++ outR where
                                        outL = (show c') ++ " : " ++ (format $ logLikelihood a) ++ " -> " ++ (format $ logLikelihood a')
                                        format = printf "%.5f"
                                        outR = replicate (len - (length outL)) ' '


getSub model a = do let piF = fromList $ safeScaledAAFrequencies a
                    case model of 
                                 WAGF -> return (wagS,piF)
                                 WAG -> return (wagS,wagPi)
                                 JTT -> return (jttS,jttPi)
                                 JTTF -> return (jttS,piF)
                                 CustomS x -> parsePamlDatIO x
                                 CustomSF x -> do (s,pi) <- parsePamlDatIO x
                                                  return (s,piF)

main = do args <- getArgs
          let ( actions, nonOpts, msgs ) = getOpt Permute options args
          opts <- foldl (>>=) (return defaultOptions) actions
          let Options {
                  optAln = aln,
                  optTree = tree,
                  optNumCats = cats,
                  optBootCount = numSim,
                  optSeed = seed,
                  optLevel = optBoot,
                  optModel = modelParams,
                  optAlg = optMethod,
                  optRaw = raw,
                  optLog = log,
                  optDebug = debugging,
                  optSub = subModel
          } = opts
          let (Logger logger) = log
          logger $ show opts
          logger $ show seed
          logger $ show numSim
          logger $ show cats 
          stdGen <- case seed of
                    Nothing -> getStdGen
                    Just x -> return $ mkStdGen x
          (aln',tree') <- case (aln,tree) of 
                                (Just a,Just t) -> do a' <- parseAlignmentFile parseUniversal a
                                                      t' <- (liftM readBiNewickTree) (readFile t)
                                                      return (a',t')
                                (_,_)           -> do printHelp opts
                                                      exitSuccess
          output <- case (aln',tree',nonOpts) of 
                (Nothing,_,_) -> do putStrLn "Failed to parse alignment"
                                    return Nothing
                (_,Left err,_) -> do putStrLn $ "Failed to parse tree " ++ err
                                     return Nothing
                (Just a,Right t,paramT)->   do logger "Debugging enabled"
                                               let method = case optMethod of 
                                                        (OptMethod a) -> a
                                                        _             -> bobyqa
                                               let pAln = pAlignment a
                                               hSetBuffering stdout NoBuffering
                                               (sMat,pi) <- getSub subModel a 
                                               let (modelF,initparams,qsetF,optF) = case modelParams of 
                                                                                Thmm a b c -> (thmmModel (cats+1) sMat pi,[a,b,c],(\x -> [toLists $ thmmModelQ (cats+1) sMat pi x]), optThmmModel method 1 cats pi sMat)
                                                                                Ras a -> (gammaModel cats sMat pi,[a],(\x -> map toLists (gammaModelQ cats sMat pi x)), optGammaModel method cats pi sMat )
                                               let (t2',priors,nClasses) = case modelParams of 
                                                                                Thmm _ _ _ -> (addModelFx (structDataN (cats+1) AminoAcid (pAln) t) (modelF initparams) (flatPriors (length $ qsetF initparams)),[1.0],(cats+1))
                                                                                Ras _ -> (addModelFx (structDataN 1 AminoAcid (pAln) t) (modelF initparams) (flatPriors (length $ qsetF initparams)),flatPriors cats,1)
                                               let nState = nClasses * 20
                                               let tol = case optBoot of
                                                                (FullOpt level) -> min level 1E-2
                                                                _               -> 1E-2
                                               (t2,params) <- case optMethod of 
                                                                                       OptNone                        -> return $ (t2',initparams)
                                                                                       (OptMethod _)                  -> do putStrLn "optimising model"
                                                                                                                            let ans = optF tol t2' initparams
                                                                                                                            let start = head ans
                                                                                                                            putStr $ getNiceOpt 100 $ head ans
                                                                                                                            mapM (\x -> do putChar '\r'
                                                                                                                                           printNiceOpt 100 x) (zip ans (tail ans))
                                                                                                                            let (a,b,_) = last ans
                                                                                                                            putStrLn ""
                                                                                                                            return (a,b)

                                               logger $ "Main model params " ++ (show params)
                                               let aX = map (fst . leftSplit) $ getAllF t2
                                               let bX = getLeftSplit t2
                                               --print $ "OK? " ++ (show (aX==bX))
                                               if (aX/=bX)
                                                       then error "Bug in smap, please report"
                                                       else return $ ()
                                               let nProc = length priors
                                               let stochmapTT lM tree = discrep $ stochmapOrder (stochmapT nProc nState (getAln tree) lM tree) (treeToMap tree) (getPriors tree)
                                               let stochmapTTA lM tree a = discrep $ stochmapOrder (stochmapT nProc nState a lM tree) (treeToMap tree) (getPriors tree)
                                               let numQuantile = 500
                                               let stdGens = take numSim $ genList stdGen
                                               let alnLength = length $ Phylo.Likelihood.columns pAln
                                               let simulations s = case (optBoot,optMethod) of 
                                                                      (FullOpt level,_) -> map (\(a,b,c) -> a) $ map (\x->last $ optF level x params) $ simulate s
                                                                      (BranchOpt,_) -> map optBLDFull0 $ simulate s
                                                                      (QuickBranchOpt,_) -> map optBLDFull0 $ simulate s
                                                                      (NoOpt,_) -> simulate s 
                                                                      where simulate x = map (\x-> patternSimulation t t2 (modelF params) priors x AminoAcid nClasses alnLength) x
                                               let simS = map (dualStochMap stochmapTT (interLMat nClasses 20) (intraLMat nClasses 20)) (simulations stdGens)
                                               let numThreads = Sync.numCapabilities
                                               let (simStochInter,simStochIntra) = if (nClasses==1)
                                                       then if (numThreads > 1)
                                                                then (undefined,(map snd simS) `using` parBuffer numThreads rdeepseq)
                                                                else (undefined,map snd simS)
                                                       else if (numThreads >1) 
                                                                then unzip (simS `using` parBuffer numThreads rdeepseq)
                                                                else unzip simS
                                               outputProgressInit 100 (fromIntegral numSim)
                                               outputProgress (const $ return ()) 100 (fromIntegral numSim) 1 simStochIntra
                                               let outputMat (name,simStochDesc,ansDesc) = do let tot x = foldr (+) (0.0) x
                                                                                              when raw $ do writeRaw (name ++"-boot-raw.txt") $ (concat . concat) simStochDesc
                                                                                                            writeRaw (name ++"-boot-edge.txt") $ concat (map (map tot) simStochDesc)
                                                                                                            writeRaw (name ++"-boot-site.txt") $ concat (map (map tot) $ map transpose simStochDesc)
                                                                                                            writeRaw (name ++"-real-raw.txt") $ concat ansDesc
                                                                                                            writeRaw (name ++"-real-edge.txt") $ map tot ansDesc
                                                                                                            writeRaw (name ++"-real-site.txt") $ map tot $ transpose ansDesc
                                                                                              let (line,lower,upper,pval) = makeQQLine numQuantile (map concat simStochDesc) (concat ansDesc)
                                                                                              let (line1,lower1,upper1,pval1) = makeQQLine numQuantile (map (map tot) simStochDesc) (map tot ansDesc)
                                                                                              let (line2,lower2,upper2,pval2) = makeQQLine numQuantile (map (map tot)  $ map transpose simStochDesc) (map tot $ transpose ansDesc)
                                                                                              let edgeQuantileMap = zip (map (\(a,b,c,d,e) -> d) $ getPartialBranchEnds t2) (perLocationQuantile numQuantile (map (map tot) simStochDesc) (map tot ansDesc))
                                                                                              let tempTree = annotateTreeWith edgeQuantileMap t2
                                                                                              writeFile (name ++ "-colour-tree.xml")  $ unlines $ quantilePhyloXML (annotateTreeWith edgeQuantileMap t2)
                                                                                              writeFile (name ++ "-pvals.txt")  $ unlines $ map (\(x,y) -> x ++ " " ++ (show y)) $ zip ["all","edge","site"] [pval,pval1,pval2]
                                                                                              m1 <- newEmptyMVar
                                                                                              forkIO $ renderableToPDFFileM m1 (makePlot line (zip lower upper) PDF) 480 480 $ name ++ "-out-all.pdf"
                                                                                              forkIO $ renderableToPDFFileM m1 (makePlot line1 (zip lower1 upper1) PDF) 480 480 $ name ++ "-out-edge.pdf"
                                                                                              forkIO $ renderableToPDFFileM m1 (makePlot line2 (zip lower2 upper2) PDF) 480 480 $ name ++ "-out-site.pdf"
                                                                                              return m1
                                               let ansIntra = stochmapTTA (intraLMat nClasses 20) t2 a
                                               let prog mVar y x = do takeMVar mVar
                                                                      putChar '\r'
                                                                      putStr $ mkProgressBar (msg "plotting") exact 100 x y
                                               m1<-outputMat ("subs",simStochIntra,ansIntra)
                                               if (nClasses==1) 
                                                       then do putStr $ mkProgressBar (msg "plotting") exact 100 0 3
                                                               mapM_ (prog m1 3) [1..3]
                                                       else do putStr $ mkProgressBar (msg "plotting") exact 100 0 6
                                                               m2<-outputMat ("switch",simStochInter,stochmapTT (interLMat nClasses 20) t2)
                                                               mapM_ (prog m1 6) [1..3]
                                                               mapM_ (prog m2 6) [4..6] 
                                               putStrLn ""
                                               putStrLn " done"
                                               return Nothing
          case output of
               Just str -> putStrLn str
               Nothing -> return ()


--outputProgress :: Int -> [a] -> IO ()
outputProgressInit width numSim = putStr $ mkProgressBar (msg "bootstrapping")  exact width 0 numSim
outputProgress f width numSim i [] = putStrLn "" --progressBar (msg "bootstrapping") exact width i numSim
outputProgress f width numSim i (x:xs) = do let m = x `seq` (mkProgressBar (msg "bootstrapping") exact width i numSim)
                                            f x 
                                            putChar '\r'
                                            putStr m
                                            outputProgress f width numSim (i+1) xs

writeRaw filename list = writeFile filename (intercalate "\n" $ map show list)

renderableToPDFFileM mvar a b c d = do renderableToPDFFile a b c d
                                       putMVar mvar ()

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
                                                fixProc pr x = foldr (+) (0.0) (map (\(x,y) -> x*y) $ zip pr x)
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
                                                                                                                                                                                                                                              
optGammaModel method cats pi s cutoff = optBSParamsBL cutoff method (0,0) (\x->0) [0.1] [Just 0.01] [Just 100.0] (flatPriors cats) model where                                                           
        model = gammaModel cats s pi
                                
optThmmModel method numModels cats pi s cutoff t2 params = optBSParamsBL cutoff method (0,numModels) (makeMapping (allIn []) t2) (map (\x->0.01) lower) lower upper [1.0] myModel t2 params where
        lower = [Just 0.01,Just 0.1,Just 0.01]                                                                                                           
        upper = [Just 0.99,Just 200.0,Just 500.0]                                                                                                                  
        myModel = thmmModel (cats+1) s pi

--this is a bit of a hack really
--all the data needed is in modelTree
--but this hack will 'recompress' the data
--to site patterns
patternSimulation origTree modelTree model priors stdGen dataType classes len = ans where
                             ans = addModelFx (structDataN classes AminoAcid (pAln) origTree) model priors 
                             pAln = pAlignment $ lAln
                             lAln = getAln $ makeSimulatedTree AminoAcid classes stdGen len modelTree 

--patternSimulation origTree modelTree model priors stdGen dataType classes len = ans where
--                             ans = makeSimulatedTree AminoAcid classes stdGen len modelTree 


getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts) where

cramerVonMises empQ theQs = pvalue where
        v = fromIntegral $ length empQ
        cram x = foldr (+) 0.0 $ map (\(w,i)->  (square $ w - (fromIntegral ((2*i)-1)/(2*v))) + (1.0/(12*v))) $ zip x (map fromIntegral [1..])
        wv2_emp = cram empQ
        wv2_thes = map cram theQs
        square x = x*x
        pvalue = (fromIntegral $ length (filter (<=wv2_emp) wv2_thes)) / (fromIntegral $ length theQs)

stochmapT nProc nState aln' lM tree = calculateStochmap nSite nState nBranch nProc nCols lM multiplicites sitemap' partials' qset' sitelikes' pi_i' branchLengths' mixProbs' where
        sitemap' = mapBack pAln'
        pAln' = pAlignment aln'
        (nSite,nCols,multiplicites) = getAlnData pAln' --nCols >= nSite by stochmap.c definition.
        nBranch =  length $ toPBELengths pBEStr
        pBEStr = getPartialBranchEnds tree
        partials' = map (map (map (toLists . trans))) $ toPBEList pBEStr
        qset' = map toLists $ getQMat tree
        sitelikes' = rawlikelihoods tree
        pi_i' = map toList $ getPi tree
        branchLengths' = toPBELengths pBEStr
        mixProbs' = getPriors tree


