{-# LANGUAGE FlexibleContexts,BangPatterns,ScopedTypeVariables #-}
import System.Environment (getArgs)
import System.Console.GetOpt
import System.Exit
import Phylo.Alignment
import Phylo.OpenBLAS as BLAS
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
import System.Unix.Directory (withTemporaryDirectory)
import System.FilePath (pathSeparator)
import Data.Binary
import Data.Maybe
import System.IO.Unsafe
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector.Unboxed as UVec
import qualified Text.JSON as JSON
import qualified Statistics.Function as SF
import qualified Data.Vector.Generic as G
import Data.Vector.Generic ((!))
import Statistics.Constants (m_epsilon)
import Control.Exception (assert)
import qualified Data.Vector.Algorithms.Tim as Tim
import Control.Monad.ST (runST)



data LMat = Inter | Intra deriving Show
data OptLevel = FullOpt Double | BranchOpt | QuickBranchOpt | NoOpt deriving Show
data SubModel = WAG | WAGF | JTT | JTTF | CustomS String | CustomSF String | JC | GTR | GTRF | HKY | HKYF deriving (Show, Eq)
whatDataType x = case x of 
                x | x `elem` nucModels -> Nucleotide
                  | otherwise  -> AminoAcid where
                                        nucModels = [JC,GTR,GTRF,HKY,HKYF]
                

--handle options like http://leiffrenzel.de/papers/commandline-options-in-haskell.html
options = [ Option ['a'] ["alignment"] (ReqArg optAlnR "FILE") "Alignment"
          , Option ['t'] ["tree"] (ReqArg optTreeR "FILE") "Tree"
          , Option ['j'] ["job-name"] (ReqArg (\arg opt -> return opt {optJobName=Just arg}) "JOBID") "set job name (required)"
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
          , Option [] ["simulate"] (ReqArg optSimulateR "LENGTH") "simulate a sequence alignment from given model"
          , Option ['p'] ["threads"] (ReqArg setThreads "THREADS") "set number of threads to use in calculation (default 1)"
            ]

setThreads arg opt = return opt {optThreads = (read arg)} 

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
optSimulateR arg opt = return opt {optSimulateOnly = Just (read arg)}

optRawR opt = return opt {optRaw = True}
optLogR opt = return opt {optLog = Logger $ hPutStrLn stderr,
                          optDebug = True}
optSubR arg opt = do let model = case arg of 
                                        "jtt" -> JTT
                                        "jtt,F" -> JTTF
                                        "wag" -> WAG
                                        "wag,F" -> WAGF
                                        "jc" -> JC
                                        "gtr" -> GTR
                                        "gtr,F" -> GTRF
                                        "hky" -> HKY
                                        "hky,F" -> HKYF
                                        x | (isF x) -> CustomSF (reverse $ drop 2 $ reverse $ x)
                                          | otherwise -> CustomS x
                     return opt {optSub = model}  where
                     isF (x:y:[]) = [x,y] == ",F"
                     isF (x:xs) = isF xs
                     isF [] = False

                                        
                                        
thmmModelR arg opt = case arg of 
                      Nothing  -> return opt {optModel = Thmm 0.1 1.0 1.0 }
                      Just args -> case (map read $ splitBy ',' args) of 
                                        [a,b,c] -> return opt { optModel = Thmm a b c }
                                        _       -> error $ "Can't parse three doubles from " ++ args

gammaModelR arg opt = case arg of 
                      Nothing -> return opt {optModel = Ras 1.0 }
                      Just args -> case (read args) of 
                                        a -> return opt { optModel = Ras a }


optNoneR opt = return opt {optAlg = OptNone}
optAlgR arg opt = case arg of 
                        Nothing -> return opt {optAlg = OptMethod bobyqa}
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
        optSub :: SubModel,
        optSimulateOnly :: Maybe Int,
        optJobName :: Maybe String,
        optThreads :: Int
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
        optBootCount = 200,
        optSeed = Nothing,
        optLevel = FullOpt 1E-1,
        optModel = Thmm 0.1 1.0 1.0,
        optAlg = OptMethod bobyqa,
        optRaw = False,
        optLog = nullOut,
        optDebug = False,
        optSub = WAGF,
        optSimulateOnly = Nothing,
        optJobName = Nothing,
        optThreads = 1
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
optFR' arg opt | (take 4 arg)=="full" = case (drop 4 arg) of 
                                                "" -> return opt {optLevel = FullOpt 1E-1}
                                                x  -> return opt {optLevel = FullOpt (read x)}
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


getSub :: SubModel -> ListAlignment -> IO SPiFunctionTuple
getSub model a = do let piF = fromList $ safeScaledAAFrequencies a
                    let nucPiF = fromList $ safeScaledNucFrequencies a
                    case model of 
                                 WAGF -> return $ customF wagS piF
                                 WAG -> return wagF
                                 JTT -> return jttF
                                 JTTF -> return $ customF jttS piF
                                 JC -> return jcF
                                 GTR -> return gtrF
                                 GTRF -> return $ (fst gtrF,zeroParam nucPiF)
                                 HKY -> return hkyF
                                 HKYF -> return $ (fst hkyF,zeroParam nucPiF)
                                 CustomS x -> do (s,pi) <- parsePamlDatIO x
                                                 return $ customF s pi
                                 CustomSF x -> do (s,pi) <- parsePamlDatIO x
                                                  return $ customF s piF

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
                  optSub = subModel,
                  optSimulateOnly = simulateOnly,
                  optJobName = jobName,
                  optThreads = numThreads
          } = opts
          let dataType = whatDataType subModel
          let prefix = case jobName of
                Nothing -> error "Please set a job name with -j"
                Just a -> a
          let (Logger logger) = prefix `seq` log
          logger $ show opts
          logger $ show seed
          logger $ show numSim
          logger $ show cats 
          stdGen <- case seed of
                    Nothing -> getStdGen
                    Just x -> return $ mkStdGen x
          tree' <- case tree of 
                       Just t          -> do t' <- (liftM readBiNewickTree) (readFile t)
                                             case t' of 
                                                Left err -> do putStrLn $ "Failed to parse tree " ++ err
                                                               exitFailure
                                                Right tx -> return tx
                       Nothing         -> do printHelp opts
                                             exitSuccess
          aln' <- case aln of 
                       Nothing -> return $ Nothing
                       Just a  -> do ans <- parseAlignmentFile parseUniversal a 
                                     case ans of
                                       Nothing -> do putStrLn "Failed to parse alignment"
                                                     exitFailure
                                       x       -> return x
          let t = tree'
          let a = case aln' of 
                        Nothing -> trace "\nERROR: You need to specify an alignment" undefined
                        Just x  -> x 
          logger "Debugging enabled"
          let debugtrace = case debugging of 
                                     True  -> trace 
                                     False -> flip const
          let method = case optMethod of 
                   (OptMethod a) -> a
                   _             -> var2
          let pAln = pAlignment a
          hSetBuffering stdout NoBuffering
          (sMatAll,piAll) <- getSub subModel a 
          let ((sMatF,sMatP,sMatL,sMatU),(piF,piP,piL,piU)) = (sMatAll,piAll)
          let extraLowerParams = sMatL++piL
          let extraUpperParams = sMatU++piU
          let sMat = (sMatF,sMatP)
          let pi = (piF,piP)
          let extraParams :: [Double] = ((getSensibleParams sMatAll)  ++ (getSensibleParams piAll))
          let alphabetSize = dataSize dataType
          let (modelF,initparams,qsetF,optF) = case modelParams of 
                                           Ras a -> (gammaModel cats sMat pi,(a:extraParams),(\x -> map toLists (gammaModelQ cats sMat pi x)), optGammaModel method cats pi sMat extraLowerParams extraUpperParams)
                                           Thmm a b c -> (thmmModel (cats+1) sMat pi,(a:b:c:extraParams),(\x -> [toLists $ thmmModelQ (cats+1) sMat pi x]), optThmmModel method 1 cats pi sMat extraLowerParams extraUpperParams)
          let (t2',priors,nClasses) = case modelParams of 
                                           Thmm _ _ _ -> (addModelFx (structDataN (cats+1) dataType (pAln) t) (modelF initparams) (flatPriors (length $ qsetF initparams)),[1.0],(cats+1))
                                           Ras _ -> (addModelFx (structDataN 1 dataType (pAln) t) (modelF initparams) (flatPriors (length $ qsetF initparams)),flatPriors cats,1)
          case simulateOnly of 
                        Just alnLength -> do let emptyAlignment = quickListAlignment (Phylo.Tree.names t) [] 
                                             let emptyTree= case modelParams of 
                                                                   Thmm a b c -> addModelFx (structDataN (cats+1) dataType (pAlignment emptyAlignment) t) (modelF initparams) (flatPriors (length $ qsetF initparams))
                                                                   Ras a      -> addModelFx (structDataN 1 dataType (pAlignment emptyAlignment) t) (modelF initparams) (flatPriors (length $ qsetF initparams))
                                             mapM putStr $ toFasta $ makeSimulatedAlignment dataType stdGen (cachedBranchModelTree emptyTree) alnLength 
                                             exitSuccess
                        Nothing        -> return ()
          let nState = nClasses * alphabetSize
          let tol = case optBoot of
                           (FullOpt level) -> min level 1E-2
                           _               -> 1E-2

          --manual control of threading
          logger $ "Set threads " ++ "1/" ++ (show numThreads) ++ "/1"

          (t2,params) <- case optMethod of 
                          OptNone         -> return $ (cachedBranchModelTree t2',initparams)
                          (OptMethod _)   -> do putStrLn "optimising model"
                                                let ans = optF tol t2' initparams
                                                let start = head ans
                                                putStr $ getNiceOpt 100 $ head ans
                                                mapM (\x -> do putChar '\r'
                                                               printNiceOpt 100 x) (zip ans (tail ans))
                                                let (a,b,_) = last ans
                                                putStrLn ""
                                                return (cachedBranchModelTree a,b)
         
          --manual control of threading
         -- Sync.setNumCapabilities (numThreads)
          let jobThreads=numThreads
          case piP of 
                0 -> logger $ "Main model params " ++ (show params)
                x -> do logger $ "Main model params " ++ (show $ reverse $ drop x $ reverse params)
                        logger $ "Nuc eqm freqs " ++ (show $ toList $ piByLog $ reverse $ take x $ reverse params)
          logger $ show t2
          let aX = map (fst . leftSplit) $ getAllF t2
          let bX = getLeftSplit t2

          --print $ "OK? " ++ (show (aX==bX))
          if (aX/=bX)
                  then error "Bug in smap, please report"
                  else return $ ()
          let nProc = length priors
          --stochmapHandle <- openFile "stochmap.out" WriteMode
          --let stochmapTTA lM tree a = discrep $ stochmapOrder (stochmapT (Just stochmapHandle) nProc nState a lM tree) (mapBack a) (getPriors tree)
          let stochmapTTA lM tree a = discrep $ stochmapOrder (stochmapT Nothing nProc nState a lM tree) (mapBack a) (getPriors tree)
          let stdGens = take numSim $ genList stdGen
          let alnLength = length $ Phylo.Likelihood.columns pAln
          let (simTrees,simAlignments) = unzip $ ((map (\x-> newSimulation a t t2 (modelF params) priors x dataType nClasses alnLength) stdGens) `using` parBuffer 3 rdeepseq)
          let simulations' = case (optBoot,optMethod) of 
                                 (FullOpt level,_) -> (map (\(a,b,c) -> debugtrace ("Params " ++ (show b)) a) $ map (\x->last $ optF level x params) $ trees,alignments)
                                 (BranchOpt,_) -> (map optBLDFull0 trees,alignments)
                                 (QuickBranchOpt,_) -> (map optBLDFull0 trees,alignments)
                                 (NoOpt,_) -> (trees,alignments)
                                 where (trees,alignments) = (simTrees,simAlignments)
          let simulations = (uncurry zip $ simulations') `using` parBuffer 1 rdeepseq
          let simS = map (dualStochMap stochmapTTA (interLMat nClasses alphabetSize) (intraLMat nClasses alphabetSize)) simulations
          let finishcalc tempdir = do 
                logger tempdir
                outputProgressInit 100 (fromIntegral numSim)
                let simS' = case jobThreads of 
                                1 -> simS
                                _ -> simS `using` parBuffer jobThreads rdeepseq
                mapM (\(i,j) -> j `deepseq` outputProgress 100 (fromIntegral numSim) i) $ zip [1..] simS'
                putStrLn "" --finish output bar
                let readIn name = do mydata <- BS.readFile name
                                     let ans = (decode mydata) :: ([[Double]],[[Double]])
                                     return ans
                let outputMat (name,fx,m1,ansDesc) = do
                                 let tot x = foldr (+) (0.0) x
                                 let d = map fx simS'

                                 when raw $ do 
                                     writeFile (name ++ "-real-raw.txt") $ annotatedSplits (map (\(a,b,c,d,e) -> d) $ getPartialBranchEnds t2) ansDesc
                                     let writeFiles x = writeFile (name ++ "-boot.raw") $ unlines $ map (annotatedSplits (map (\(a,b,c,d,e) -> d) $ getPartialBranchEnds t2)) x
                                     withFile (name ++ "-boot-raw.txt") WriteMode $ (\fh ->
                                        do let mydata = d 
                                           let mydata'  = map (annotatedSplits (map (\(a,b,c,d,e) -> d) $ getPartialBranchEnds t2)) mydata
                                           mapM (hPutStr fh) mydata' )
                                     return ()
                                 let qqFunc f2 proc x = do
                                                      let boots = map f2 x
                                                      let real = f2 ansDesc
                                                      when (raw && proc/="raw") $ do
                                                                    writeRaw (name ++"-boot-" ++ proc ++".txt") $ concat boots
                                                                    writeRaw (name ++"-real-" ++ proc ++".txt") $ real
                                                      return $ makeQQLine boots real
                                 (line,lower,upper,pval) <- qqFunc concat "raw"  d
                                 renderableToPDFFile (makePlot line (zip lower upper) PDF) 480 480 $ name ++ "-all.pdf"
                                 putMVar m1 ()
                                 (line1,lower1,upper1,pval1) <- qqFunc (map tot) "edge" d
                                 renderableToPDFFile (makePlot line1 (zip lower1 upper1) PDF) 480 480 $ name ++ "-edge.pdf"
                                 putMVar m1 ()
                                 (line2,lower2,upper2,pval2) <- qqFunc (map tot . transpose) "site" d
                                 renderableToPDFFile (makePlot line2 (zip lower2 upper2) PDF) 480 480 $ name ++ "-site.pdf"
                                 putMVar m1 ()
                                 let edgeQuantileMap x = zip (map (\(a,b,c,d,e) -> d) $ getPartialBranchEnds t2) (perLocationQuantile (map (map tot) x) (map tot ansDesc))
                                 let eQM = edgeQuantileMap d
                                 let tempTree = annotateTreeWith eQM t2
                                 writeFile (name ++ "-colour-tree.xml")  $ unlines $ quantilePhyloXML (annotateTreeWith eQM t2)
                                 putMVar m1 ()
                                 writeFile (name ++ "-pvals.txt")  $ unlines $ map (\(x,y) -> x ++ " " ++ (show y)) $ zip ["all","edge","site"] [pval,pval1,pval2]
                                 putMVar m1 ()
                let ansIntra = stochmapTTA (intraLMat nClasses alphabetSize) t2 pAln
                m1<-newEmptyMVar
                forkIO $ outputMat (prefix++"-subs",snd,m1,ansIntra)
                let prog mVar y x = do 
                    putChar '\r'
                    takeMVar mVar
                    putStr $ mkProgressBar (msg "plotting    ") exact 100 x y
                if nClasses==1
                   then do putStr $ mkProgressBar (msg "plotting    ") exact 100 0 5
                           mapM_ (prog m1 5) [1..5]
                   else do putStr $ mkProgressBar (msg "plotting    ") exact 100 0 10
                           forkIO $ outputMat (prefix++"-switch",fst,m1,stochmapTTA (interLMat nClasses alphabetSize) t2 pAln)
                           mapM_ (prog m1 10) [1..10]
                putStrLn ""
                putStrLn " done"
                return Nothing
          withTemporaryDirectory "smap" finishcalc 


annotatedSplits splits columns = unlines $ map s $ zip (map annotate splits) (map p columns) where
        annotate (a,b) = (showList "," a) ++ " " ++ (showList "," b)
        showList x = foldr (++) "" . intersperse x 
        p = showList " " . map show
        s (a,b) = a ++ " " ++ b


--outputProgress :: Int -> [a] -> IO ()
outputProgressInit width numSim = putStr $ mkProgressBar (msg "bootstrapping")  exact width 0 numSim
outputProgress width numSim i = do 
   putChar '\r'
   putStr $ mkProgressBar (msg "bootstrapping") exact width i numSim

writeRaw filename list = writeFile filename (intercalate "\n" $ map show list)


dualStochMap stochmapTT lM1 lM2 (t,a) = (stochmapTT lM1 t a, stochmapTT lM2 t a) 
length2 = map length
length3 = map length2

--FIXME not sure about this
reverseQuantile myQList point = case (filter (\(x,y) -> point <= y) $ zip [0..] myQList) of
                                        (x:xs) -> fst $x
                                        [] -> length $ myQList

discrep :: ([[Double]], [Double]) -> [[Double]]
discrep ((x:xs),(y:ys)) = (map (\i->i-y) x):(discrep (xs,ys))
discrep ([],[]) = []

perLocationQuantile :: [[Double]] -> [Double] -> [Double]
perLocationQuantile simulated locdist = map ((/ (fromIntegral numQuantile)) . fromIntegral . (reverseQuantile simQuantiles)) locdist where
        numQuantile = length locdist
        quantileF d q = continuousBy medianUnbiased q numQuantile d
        simVec = UVec.fromList $ concat simulated
        simQuantiles = map (quantileF simVec) [0..numQuantile]



-- number of Quantiles needed
-- simulated data [sim_n [edge/site]]

makeQQLine :: [[Double]] -> [Double] -> ([Double],[Double],[Double],Double) 
makeQQLine simulated empirical = (empQuantile,lower,upper,pvalue) where
        simcount = length simulated
        numQ = simcount
        f (real,sim) = (fromIntegral $ length $ filter (<= real) sim) / (fromIntegral simcount)

        fval x = timSort $ map f $ zip x (transpose simulated)
        fvals = fval empirical
        fvalSim = map fval simulated

        makeAllQ = continuousAll medianUnbiased numQ 
        empQuantile = makeAllQ fvals 
        simQuantile = map makeAllQ fvalSim
        simQuantile' = map timSort $ transpose simQuantile
        
        pvalue = 1.0 - (cramerVonMises empQuantile simQuantile)
        lowerI = floor $ (fromIntegral simcount) * 0.025
        upperI = ceiling ((fromIntegral simcount) * 0.975) -1
        lower = map (UVec.! lowerI) simQuantile'
        upper = map (UVec.! upperI) simQuantile'

timSort list = runST $ do 
  vec <- UVec.thaw $ UVec.fromList list
  Tim.sort vec
  UVec.freeze vec


continuousAll' cParam q x = map (\k -> continuousBy cParam k q x) [0..q]

{-# INLINE continuousAll #-}
continuousAll :: G.Vector v Double =>
                 ContParam  -- ^ Parameters /a/ and /b/.
              -> Int        -- ^ /q/, the number of quantiles.
              -> v Double   -- ^ /x/, the sample data.
              -> [Double]
continuousAll (ContParam a b) q x = 
    assert (q >= 2) .
    assert (G.all (not . isNaN) x) $
    map pFunc [0..q]
  where
    sx              = x
    bracket m       = min (max m 0) (n - 1)
    eps             = m_epsilon * 4
    n               = G.length x
    pFunc k = (1-h) * item (j-1) + h * item j where 
      j               = floor (t + eps)
      t               = a + p * (fromIntegral n + 1 - a - b)
      p               = fromIntegral k / fromIntegral q
      h | abs r < eps = 0
        | otherwise   = r
        where r       = t - fromIntegral j
      item            = (sx !) . bracket


--makeQQLine :: Int -> [[Double]] -> [Double] -> ([Double],[Double],[Double],Double)
--makeQQLine numQuantile simulated empirical = (empQuantile,lower,upper,pvalue) where
--                                        revQuantile x = map (quantileF (revQuantileRawF x)) [0..numQuantile]
--                                    
--                                        revQuantileRawF x = UVec.fromList $ map (\y->(fromIntegral $ reverseQuantile simQuantiles y)/(fromIntegral numQuantile)) x
--
--                                        quantileF d q = continuousBy medianUnbiased q numQuantile d
--                                        tSim = transpose simulated
--                                        simVec i = UVec.fromList $ (tSim !! i)
--                                        simQuantiles = map (\x -> quantileF (simVec x) x) [0..numQuantile]
--
--                                        numSim = length simulated
--                                        revQuantileSimS' = map revQuantile simulated
--                                        revQuantileSimS = map sort $ transpose revQuantileSimS'
--                                        empQuantile = revQuantile empirical
--                                        pvalue = 1.0 - (cramerVonMises empQuantile revQuantileSimS')
--
--                                        lowerI = floor $ (fromIntegral numSim) * 0.025
--                                        upperI = (ceiling ((fromIntegral numSim) * 0.975)) - 1 
--                                        lower = map (!!lowerI) revQuantileSimS
--                                        upper = map (!!upperI) revQuantileSimS
--


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

safeScaledFrequencies f = normalise . map (\x-> if x < 1e-15 then 1e-15 else x) . f
safeScaledAAFrequencies = safeScaledFrequencies scaledAAFrequencies
safeScaledNucFrequencies = safeScaledFrequencies $ scaledFrequencies "TCAG"

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
                                                fixProc pr x = foldl' (+) (0.0) (map (\(x,y) -> x*y) $ zip pr x)
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
                                                 total xs = foldl' (+) 0.0 xs
                                                 fmtBranch2 bs [] = []
                                                 fmtBranch2 (b:bs) ((cond,prior):xs) = ((show b) ++"\t" ++ (show $ total cond) ++ "\t" ++ (show $ total prior) ++ "\n") : (fmtBranch2 bs xs)
                                                 headerSite = "Site\tTotal_conditional_expectation\tTotal_prior_expectation\n"
                                                 remainderSite = fmtBranch2 [0..] $ zip condE'' priorE''

qList qFunc numQ = map (\x-> qFunc x numQ) [0..numQ]




--uniformQ :: [([[Double]],[Double])] -> (Int -> Int -> Double)
--uniformQ simData = partialFunc where 
--                      partialFunc q k | trace ("Made vec " ++ (show vec)) True = continuousBy medianUnbiased q k vec
--                      vec | trace ("Making vector " ++ (take 10 $ show simData)) True = (UVec.fromList $ concatMap allDisc simData)

reformatAns mapping priors stochResult = (map $ uncurry3 stochmapOrder) $ zip3 stochResult mapping priors
--uniformQRaw mapping priors stochResult = uniformQ $ reformatAns mapping priors stochResult

uncurry3 f (a,b,c) = f a b c 

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
                                                                                                                                                                                                                                              
optGammaModel method cats pi s lower upper cutoff = optBSParamsBL cutoff method (0,0) (\x->0) (replicate ((length lower)+1) 0.1) ((Just 0.01):lower) ((Just 100.0):upper) (flatPriors cats) model where                                                           
        model = gammaModel cats s pi
                                
optThmmModel method numModels cats pi s lower' upper' cutoff t2 params = optBSParamsBL cutoff method (0,numModels) (makeMapping (allIn []) t2) (map (\x->0.01) lower) lower upper [1.0] myModel t2 params where
        lower = [Just 0.01,Just 0.1,Just 0.01] ++ lower'                                                                                                       
        upper = [Just 0.99,Just 200.0,Just 500.0] ++ upper'                                                                                                           
        myModel = thmmModel (cats+1) s pi

newSimulation origAln origTree tree model priors stdGen dataType classes len = (addModelFx (structDataN classes dataType (pAln) origTree) model priors,pAln) where
        aln = makeSimulatedAlignmentWithGaps dataType stdGen tree origAln
        pAln = pAlignment aln


getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts) 

cramerVonMises empQ theQs = pvalue where
        v = fromIntegral $ length empQ
        cram x = foldl' (+) 0.0 $ map (\(w,i)->  (square $ w - (fromIntegral ((2*i)-1)/(2*v))) + (1.0/(12*v))) $ zip x (map fromIntegral [1..])
        wv2_emp = cram empQ
        wv2_thes = map cram theQs
        square x = x*x
        pvalue = (fromIntegral $ length (filter (<=wv2_emp) wv2_thes)) / (fromIntegral $ length theQs)

stochmapT fh nProc nState pAln' lM tree = unsafePerformIO $ calculateAndWrite nSite nState nBranch nProc nCols lM multiplicites sitemap' partials' qset' sitelikes' pi_i' branchLengths' mixProbs' fh where
        sitemap' = mapBack pAln'
        (nSite,nCols,multiplicites) = getAlnData pAln' --nCols >= nSite by stochmap.c definition.
        nBranch =  length $ toPBELengths pBEStr
        pBEStr = getPartialBranchEnds tree
        partials' = map (map (map (toLists . trans))) $ toPBEList pBEStr
        qset' = map toLists $ getQMat tree
        sitelikes' = rawlikelihoods tree
        pi_i' = map toList $ getPi tree
        branchLengths' = toPBELengths pBEStr
        mixProbs' = getPriors tree


