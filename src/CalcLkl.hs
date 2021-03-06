import System.Environment (getArgs)
import System.Console.GetOpt
import Phylo.Alignment
import Data.List
import Control.Monad
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
import Phylo.NeXML
import Phylo.NLOpt
import Phylo.OpenBLAS as BLAS
import Phylo.OptOutput
import System.Exit (exitSuccess,exitFailure)

data Flag = JobName String | PerBranch String | NumCats String | AlignmentFile String | TreeFile String | Alpha String | Sigma String | Opt | OptAlpha | PriorZero String | OptThmm | OptThmmP | OptThmm2 String | OptSim1 | OptSim2 String | OptSim0 String | Seed String | ThmmStochmap String | NoOpt 
        deriving (Show,Eq,Ord)
options = [ Option ['a'] ["alignment"] (ReqArg (\arg opt -> return opt {optAln=Just arg})   "FILE") "Alignment",
            Option ['t'] ["tree"] (ReqArg (\a o -> return o {optTree=Just a}) "FILE") "Tree",
            Option ['j'] ["job-name"] (ReqArg (\a o -> return o {optJobName=Just a}) "PREFIX") "Job Name",
            Option ['s'] ["seed"] (ReqArg (\a o -> return o {optSeed = Just $ read a})  "INTEGER") "Seed",
            Option [] ["num-cats"] (ReqArg (\a o -> return o {optNumCats=read a}) "NUMGAMMACATS") "Number of Gamma Rate Categories",
            Option ['g'] ["gamma"] (ReqArg (\a o -> return o {optAlpha=read a}) "DOUBLE") "Use Gamma Model" ,
            Option ['S'] ["sigma"] (ReqArg (\a o -> return o {optSigma=Just $ read a}) "DOUBLE") "Use Switching THMM", 
            Option ['z'] ["prior-zero"] (ReqArg (\a o -> return o {optPriorZero=Just $ read a}) "DOUBLE") "Set Prior for zero rate in THMM",
            Option ['o'] ["opt"] (NoArg (\o -> return o {optOpt=True})) "Optimise Model",
            Option [] ["per-branch"] (ReqArg (\a o -> return o {optPerBranch=Just a}) "splits") "Per-branch sigma",
            Option [] ["all-branch"] (NoArg (\o -> return o {optAllBranch=True})) "Sigma for each branch"]


data Options = Options  {
        optAln  :: Maybe String,
        optTree :: Maybe String,
        optJobName :: Maybe String,
        optNumCats :: Int,
        optSeed :: Maybe Int,
        optOpt :: Bool,
        optSigma :: Maybe Double,
        optPriorZero :: Maybe Double,
        optAlpha :: Double,
        optPerBranch :: Maybe String,
        optAllBranch :: Bool
} deriving Show

defaultOptions = Options{
        optAln=Nothing,
        optTree=Nothing,
        optJobName=Nothing,
        optNumCats=4,
        optSeed=Nothing,
        optOpt=False,
        optSigma=Nothing,
        optPriorZero=Nothing,
        optAlpha=1.0,
        optPerBranch=Nothing,
        optAllBranch=False
}

wagSX = ((\_->wagS),0)

wagPiF = piX wagPi
piX piF = ((\_->piF),0)

getMapping t spl = ((1,(length goodNodes)+1),Just $ makeMapping (allIn goodNodes) t) where
             goodNodes = map (splitBy ',') $ splitBy ' ' spl

printHelp' exit opt = do putStrLn $ usageInfo header options
                         putStrLn example
                         exit where
                            header = unlines ["AmpliPhy: Optimise THMM and Gamma models"]
                            example = unlines []

printHelp = printHelp' exitSuccess
printHelpError = printHelp' exitFailure


main = do args <- getArgs
          BLAS.set_num_threads 1
          let ( actions, nonOpts, msgs ) = getOpt Permute options args
          opts <- foldl (>>=) (return defaultOptions) actions   
          let helpOutput="HELP"
          let Options {
                  optAln=alnMaybe,
                  optTree=treeMaybe,
                  optJobName=jobNameMaybe,
                  optNumCats=cats,
                  optSeed=seedMaybe,
                  optOpt=optFlag,
                  optSigma=sigmaMaybe,
                  optPriorZero=priorZero'',
                  optAlpha=alpha',
                  optPerBranch=mappingMaybe,
                  optAllBranch=optAllSwitch
          } = opts
          -- we always do the zero rate matrix in a THMM
          -- but with no switching, if the user specifies -z
          -- we should include a zero rate matrix
          -- otherwise, don't (just do a gamma model)
          let (priorZero',doZeroGamma) = case priorZero'' of  
                Nothing -> (0.05,False)
                Just x  -> (x,True)
          aln <- case alnMaybe of
                Just a  -> parseAlignmentFile parseUniversal a
                Nothing -> do putStrLn "Please specify an alignment"
                              printHelpError opts
          tree <- case treeMaybe of
                Just a  -> (liftM readBiNewickTree) (readFile a)
                Nothing -> do putStrLn "Please specify a tree" 
                              printHelpError opts
          stdGen <- case seedMaybe of 
                Just a  -> return $ mkStdGen a
                Nothing -> getStdGen
          jobName <- case jobNameMaybe of 
                Just a  -> return a
                Nothing -> do putStrLn "Please specify a job name"
                              printHelpError opts
          output <- case (aln,tree,sigmaMaybe) of 
                (Just a,Right t,Nothing) ->               do let shouldOpt = optFlag
                                                             (optTree,params) <- case shouldOpt of 
                                                                                       True ->   optWithOutput progress
                                                                                       False ->  return (addModelNNode t2 m' pri pi',initParams) where
                                                                                                     (m',pi',pri) = model initParams
                                                             let [priorZero,alpha] = case doZeroGamma of
                                                                        False -> [0.0,(head params)]
                                                                        True  -> params
                                                             let lkl = logLikelihood optTree
                                                             writeFile (jobName ++ "-out.tre") (show optTree)
                                                             writeFile (jobName ++ "-out.lkl") ((show lkl) ++ "\n")
                                                             case doZeroGamma of
                                                                        False -> do writeFile (jobName ++ "-out.param") ("alpha = " ++ (show alpha) ++ "\n")
                                                                        True -> writeFile (jobName ++ "-out.param") ("alpha = " ++ (show alpha) ++ "\npriorZero = " ++ (show priorZero) ++ "\n")
                                                             return $ Just $ "Gamma: " ++ (show alpha) ++ " priorZero: " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                                 t2 = structDataN 1.0 1 AminoAcid (pAlignment a) t
                                                                 piF = fromList $ safeScaledAAFrequencies a
                                                                 model = case doZeroGamma of 
                                                                        False -> gammaModel cats wagSX $ piX piF
                                                                        True  -> zeroGammaModel cats wagSX $ piX piF 
                                                                 progress = optParamsAndBLIO model t2 initParams lBound uBound 0.01
                                                                 (initParams,lBound,uBound) = case doZeroGamma of 
                                                                     True -> ([priorZero',alpha'],[Just 0.0,Just 0.001],[Just 0.99,Nothing])
                                                                     False -> ([alpha'],[Just 0.001],[Nothing])

                (Just a,Right t,Just sigma') -> do 
                        --let scaleT = 1e300
                        let scaleT = 1
                        let (scale,priorScale) = getScales t scaleT
                        let priors = [1.0]
                        let t2 = structDataN scale (cats+1) AminoAcid (pAlignment a) t                                                                                                             
                        let piF = fromList $ normalise $ map (\x-> if x < 1e-15 then 1e-15 else x) $ safeScaledAAFrequencies a
                        let (numBSParams,mapping) = case mappingMaybe of
                                Nothing -> ((0,0),Nothing)
                                Just s  -> getMapping t2 s
                        let model = thmmModel (cats+1) wagSX (piX piF)
                        let progress = optBSParamsAndBLIO (0,0) Nothing model t2 [priorZero',alpha',sigma'] [Just 0.01,Just 0.001, Just 0.00] [Just 0.99,Just 100.0,Just 10000.0] 0.01
                        let shouldOpt = numBSParams /= (0,0) || optFlag --always opt if numBSParams /= 0
                        (optTree',[priorZero,alpha,sigma]) <- case shouldOpt of 
                                                                True -> optWithOutput progress
                                                                False -> return (addModelNNode t2 m' priors pi',[priorZero',alpha',sigma']) where
                                                                           (m',pi',_) = model [priorZero',alpha',sigma']
                        ans <- case (numBSParams,optAllSwitch) of
                             (_,True) -> do --every branch has a sigma
                                let model = thmmPerBranchModel (cats+1) wagS piF
                                let progress = trace ("Per Branch! " ) $ optBSParamsAndBLIO numBSParams mapping model (annotateTreeWith' (\x->sigma) optTree') [priorZero,alpha] [Just 0.01,Just 0.001] [Just 0.99,Just 100.0] 0.01
                                (optTree2,multiParams) <- optWithOutput progress
                                let optTree = annotateTreeWithNumberSwitches priorScale AminoAcid optTree2
                                writeFile (jobName ++ "-out.tre") $ (show $ getSubBL 0 optTree) ++ "\n"
                                writeFile (jobName ++ "-switching-out.tre") $ (show $ getSubBL 2 optTree) ++ "\n"
                                writeFile (jobName ++ "-out.lkl") ((show $logLikelihood optTree2) ++ "\n")
                                let sigmas = map (!!1) $ getBL' optTree2 []
                                writeFile (jobName ++ "-out.param") $ "alpha = " ++ (show alpha) ++ "\np_0 = " ++ (show priorZero) ++ "\nsigma = " ++ (show sigmas) ++"\n"
                                print multiParams
                                print optTree
                                return $ Just "OK"
                             ((0,0),False) -> do --just one sigma
                                let optTree = annotateTreeWithNumberSwitchesSigma priorScale AminoAcid sigma optTree'
                                let lkl = logLikelihood optTree
                                writeFile (jobName ++ "-out.tre") $ (show $ getSubBL 0 optTree) ++"\n"
                                writeFile (jobName ++ "-switching-out.tre") $ (show $ getSubBL 2 optTree) ++ "\n"
                                writeFile (jobName ++ "-out.param") $ "alpha = " ++ (show alpha) ++ "\np_0 = " ++ (show priorZero) ++ "\nsigma = " ++ (show sigma) ++ "\n"
                                writeFile (jobName ++ "-out.lkl") ((show $logLikelihood optTree') ++ "\n")
                                return $ Just $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show sigma) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) ++ "\n" ++ (show optTree) 
                             (_,False) -> do  --hybrid
                                let model = thmmPerBranchModel (cats+1) wagS piF
                                let (a,b) = numBSParams
                                let progress = trace ("Per Branch! " ++ (show [a,b])) $ optBSParamsAndBLIO numBSParams mapping model optTree' ((replicate (a*b) sigma)  ++ [priorZero,alpha]) ((replicate (a*b) $ Just 0.00) ++ [Just 0.01,Just 0.001]) ((replicate (a*b) $ Just 10000.0) ++ [Just 0.99,Just 100.0]) 0.01
                                (optTree2,multiParams) <- optWithOutput progress
                                let (alpha:priorZero:sigmas) = reverse multiParams 
                                {- sigmas is now reversed
                                   internal convention is that for --per-branch="set1 set2"
                                   sigmas = [sigma1,sigma2,sigma_background]
                                   makes most sense to output [sigma_background,sigma1,sigma2] I think --}
                                let sigmas' = ((head sigmas):(reverse $ tail sigmas))
                                let optTree = annotateTreeWithNumberSwitches priorScale AminoAcid optTree2
                                let lkl = logLikelihood optTree
                                writeFile (jobName ++ "-out.tre") $ (show $ getSubBL 0 optTree) ++ "\n"
                                writeFile (jobName ++ "-switching-out.tre") $ (show $ getSubBL 2 optTree) ++ "\n"
                                writeFile (jobName ++ "-out.param") $ "alpha = " ++ (show alpha) ++ "\np_0 = " ++ (show priorZero) ++ "\nsigma = " ++ (show sigmas') ++"\n"
                                writeFile (jobName ++ "-out.lkl") ((show $logLikelihood optTree2) ++ "\n")
                                return $ Just $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show sigmas) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) ++ "\n" ++ (show optTree) 
                        return ans
{--                (Just a,Right t,(ThmmStochmap params ):[])-> do let alpha:sigma:priorZero:[] = map read $ take 3 $ words params
                                                                let interintra = head $ drop 3 $ words params
                                                                let pAln = pAlignment a
                                                                let (nSite,nCols,multiplicites) = getAlnData pAln where
                                                                    getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts)
                                                                let piF = fromList $ safeScaledAAFrequencies a
                                                                let model = thmmModel (cats+1) wagS piF [priorZero,alpha,sigma]
                                                                let t2 = addModelFx (structDataN (cats+1) AminoAcid (pAln) t) model [1.0]
                                                                let nState = (cats+1)*20
                                                                let lMat = case interintra of 
                                                                                "inter" -> interLMat (cats+1) 20
                                                                                "intra" -> intraLMat (cats+1) 20
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
                                                                --handle <- openFile "out.txt" WriteMode
                                                                putStrLn $ show $ zip [0..] $ getCols pAln 
                                                                putStr $ splitsStr $ toPBESplits pBEStr
                                                                --putStrLn $ "OK " ++ (show $ accEigQ (qMat) ((snd model)!!0))
                                                                let xml = xmlTree t2
                                                                putStrLn "Here comes the xml:"
                                                                putStr xml
                                                                ans <- calculateAndWrite nSite nState nBranch nProc nCols lMat multiplicites sitemap partials qset sitelikes pi_i branchLengths mixProbs Nothing
                                                                print sitemap
                                                                stochmapOut ans sitemap [1.0] putStr
                                                                print "OK"
                                                                return Nothing
                                                                {-c_test chandle $ fromIntegral 10-}
                (Just a,Right t,(OptThmmP):[])-> return $ Just $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show optTree) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        piF = fromList $ safeScaledAAFrequencies a
                                                        t2 = addModelFx (structDataN (cats+1) AminoAcid (pAlignment a) t) (gammaModel cats wagS piF [0.5]) $ flatPriors cats                                                                                                                                                                
                                                        startBL = getBL t2
                                                        newBL = map (\x -> x:0.0:[]) startBL
                                                        t3 = fst $ setBLX' 0 newBL t2
                                                        model = thmmPerBranchModel (cats+1) wagS piF
                                                        (optTree,[priorZero,alpha]) = optParamsAndBL model t3 [0.1,0.5] [1.0] [Just 0.01,Just 0.001] [Just 0.99,Nothing] 0.01
                                                        bls = getBL t3
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(OptThmm2 spl'):xs)-> do (optTree',optParams) <- case xs of 
                                                                list | NoOpt `elem` list -> return (t4,p) where
                                                                                              p = (initSigma ++ [initP0,initAlpha])
                                                                                              t4 | trace ("NUMMODELS " ++ (show numModels)) True = getFuncT1A [1.0] mapped (1,numModels) model t3 p
                                                                _ -> optBSParamsBL 1E-2 bobyqa (1,numModels) mapped 4E-1
                                                         let optTree = annotateTreeWithNumberSwitches optTree'
                                                         print optTree
                                                         let (sigma,[priorZero,alpha]) = splitAt numModels optParams
                                                         let lkl = logLikelihood optTree
                                                         putStrLn ("SIGMA " ++ (show sigma))
                                                         --
                                                         let (nSite,nCols,multiplicites) = getAlnData (pAlignment a) where
                                                                    getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts)
                                                         let nState = (cats+1)*20
                                                         let pi_i = map toList $ getpi optTree' where
                                                                        getpi (DTree _ _ _ _ _ _ pi) = pi
                                                         let sitelikes = [likelihoods optTree']
                                                         let nProc = 1 {-- fixme --}
                                                         let pBEStr = getPartialBranchEnds optTree'
                                                         let qset = map toLists $ map head $ toPBEQ pBEStr
                                                         let mixProbs = getPriors optTree'
                                                         let nBranch =  length $ toPBELengths pBEStr
                                                         let branchLengths = toPBELengths pBEStr
                                                         let branchLengths2 = everyOther $ getBL optTree' where 
                                                                everyOther (x:y:z) = x:(everyOther z)
                                                                everyOther [] = []
                                                         print branchLengths 
                                                         print branchLengths2
                                                         --TODO TEST
                                                         let interintra = "inter"
                                                         let sitemap = mapBack $ pAlignment a
                                                         let partials = map (map (map (toLists . trans))) $ toPBEList pBEStr
                                                         let lMat = case interintra of 
                                                                   "inter" -> interLMat (cats+1) 20
                                                                   "intra" -> intraLMat (cats+1) 20
                                                         putStrLn $  "Thmm: " ++ (show alpha) ++ " " ++ (joinWith " " $ sigma) ++ " " ++ " " ++ (show priorZero) ++ " " ++ (show lkl) ++ "\n" ++ (show optTree) 
                                                         putStrLn $ posteriorTipsCSV (cats+1) (pAlignment a) optTree'
                                                         --calculateAndWrite nSite nState nBranch nProc nCols lMat multiplicites sitemap partials qset sitelikes pi_i branchLengths mixProbs stdout
                                                         return Nothing where
                                                                spl = clean spl'
                                                                piF = fromList $ safeScaledAAFrequencies a
                                                                t2 = addModelFx (structDataN (cats+1) AminoAcid (pAlignment a) t) (gammaModel cats wagS piF [0.5]) $ flatPriors cats                                                                                                                                                             
                                                                startBL = getBL t2
                                                                newBL = map (\x -> x:0.0:[]) startBL
                                                                t3 = fst $ setBLX' 0 newBL t2
                                                                goodNodes = map (splitBy ',') $ drop 1 $ splitBy ' ' spl
                                                                mapped = makeMapping (allIn goodNodes) t3
                                                                numModels = length $ nub $ whichModels t3 mapped
                                                                initParams = map read $ splitBy ',' $ head $ splitBy ' ' spl
                                                                initAlpha = head initParams
                                                                initSigma = take numModels $ tail initParams
                                                                initP0 = last initParams
                                                                model | traceShow t3 True = thmmPerBranchModel (cats+1) wagS piF 
                                                                lower = (replicate numModels $ Just 0.0) ++ [Just 0.001,Just 0.001] 
                                                                upper = (replicate numModels Nothing) ++ [Just 0.99,Nothing]
                                                        {-optBSParamsBL numBSParam mapping = optWithBS' [] 1E-4 numBSParam (Just mapping)-}
                                                        {-optParamsBL = optWithBS' [] 1E-4 (0,0) Nothing-}
                                                        {-optWithBS' :: [(Double,IterType)] -> Double -> (Int,Int) -> (Maybe (DNode -> Int)) -> [Maybe Double] -> [Maybe Double] -> [Double] -> ModelF -> DNode -> [Double] -> (DNode,[Double])-}
                                                        {-optWithBS' iterations cutoff numBSParam mapping lower upper priors model tree startParams = (bestTree,bestParams) where-}
{-                (Just a,Right t,OptSim1:[]) -> putStrLn $ concat $ toFasta $ simulateSequences AminoAcid (cats+1) stdGen 349 t2 where-}
                                                        {-alpha = 0.8335109218715334-}
                                                        {-priorZero = 0.1907077998691572-}
                                                        {-sigma0 = 15.631076379213784-}
                                                        {-sigma1 = 1.0836675920110694e-8-}
                                                        {-piF = fromList $ scaledAAFrequencies a-}
                                                        {-t2 = addModelFx (setBLMapped 1 (dummyTree (structDataN (cats+1) AminoAcid (pAlignment a) t)) mapped ) (thmmPerBranchModel (cats+1) cpRevS cpRevPi [priorZero,alpha]) [1.0]-}
                                                        {-goodNodes = ["E_Nosloc","E_Enccun","E_Gluple","A_Aerper","A_Metbar"] -}
                                                        {-mapped = makeMapping (\(x,y) -> if (((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))) then ([sigma0]) else ([sigma1])) t2-}
                (a,Right t,(OptSim0 params):[]) -> return $ Just $ concat $ toFasta $ makeSimulatedAlignment AminoAcid stdGen t2 (floor length)  where  
                                                        alpha:sigma0:sigma1:priorZero:length:[] = map (read) $take 5 $ words params
                                                        goodNodes = splitBy ',' $ head $ drop 5 $ words params 
                                                        piF = case a of 
                                                                Just aln -> fromList $ safeScaledAAFrequencies aln
                                                                Nothing -> wagPi
                                                     --   piF | trace ((show sigma0) ++ "  " ++  (show sigma1)) True= fromList $ scaledAAFrequencies a
                                                        emptyAlignment = ListAlignment [] [] []
                                                        t2 = trace (show goodNodes) $ addModelFx (setBLMapped 1 (dummyTree (structDataN (cats+1) AminoAcid (pAlignment emptyAlignment) t)) mapped ) (thmmPerBranchModel (cats+1) wagS piF [priorZero,alpha]) [1.0]
                                                        mapped = makeMapping (\(x,y) -> if (((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))) then ([sigma0]) else ([sigma1])) t2
                (Just a,Right t,OptSim2 params:[]) -> do (lgInnerS,lgInnerPi) <- parsePamlDatIO "lgInner"
                                                         (lgOuterS,lgOuterPi) <- parsePamlDatIO "lgOuter"
                                                         return $ Just $ compute (lgInnerS,lgInnerPi) (lgOuterS,lgOuterPi) where
                                                             compute (s0,pi0) (s1,pi1) = concat $ toFasta $ makeSimulatedAlignment AminoAcid stdGen t3 349 where
                                                                alpha:sigma0:sigma1:priorZero:[] = map (read) $ words params
                                                                piF = fromList $ safeScaledAAFrequencies a
                                                                (model0,pi_0) = thmmPerBranchModel (cats+1) s0 pi0 [priorZero,alpha]
                                                                (model1,pi_1) = thmmPerBranchModel (cats+1) s1 pi1 [priorZero,alpha]
                                                                t2 = addModelFx (setBLMapped 1 (dummyTree (structDataN (cats+1) AminoAcid (pAlignment a) t)) mapped ) (thmmPerBranchModel (cats+1) cpRevS cpRevPi [priorZero,alpha]) [1.0]
                                                                tfFunc x y = ((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))
                                                                goodNodes = ["E_Nosloc","E_Enccun","E_Gluple","A_Aerper","A_Metbar"] 
                                                                mapped = makeMapping (\(x,y) -> if (tfFunc x y)  then ([sigma0]) else ([sigma1])) t2
                                                                mappedModels = makeMapping (\(x,y) -> if (tfFunc x y) then model1 else model0) t2
                                                                t3 = restructDataMapped t2 mappedModels [1.0] pi_0
                                                                --}
                (_,_,_) -> error "Can't parse something"
          case output of
               Just str -> putStrLn str
               Nothing -> return ()


{-which of sets am I in (for split (l,r)? or (length sets) if I am not-}
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
                 total = foldl' (+) 0.0 list

safeScaledAAFrequencies = normalise . map (\x-> if x < 1e-15 then 1e-15 else x) . scaledAAFrequencies
optBSParamsAndBLIO numBSParam mapping model tree params lower upper cutoff = optWithBSIO' bobyqa [] cutoff numBSParam mapping (map (\x->0.01) lower) (map (\x->1E-4) lower) lower upper model (dummyTree tree) params 
optParamsAndBLIO = optBSParamsAndBLIO (0,0) Nothing
trim = f . f where 
   f = reverse . dropWhile isSpace
clean = clean' . trim where
   clean' (' ':' ':xs) = clean' (' ':xs)
   clean' (x:xs) = x:(clean' xs)
   clean' [] = []
joinWith i (x:x':xs) = (show x) ++ i ++ (joinWith i (x':xs))
joinWith i [x] = (show x) 
joinWith i [] = []

stochmapOrder :: ([[[Double]]],[[Double]]) -> [Int] -> [Double] -> ([[Double]],[[Double]])
stochmapOrder (condE,priorE) mapping priors = order where
                                                condE' = map (fixProc2 priors) condE
                                                priorE' = map (fixProc priors) priorE
                                                fixProc pr x = (foldl' (+) (0.0) (map (\(x,y) -> x*y) $ zip pr x))
                                                fixProc2 pr xs = (map $ fixProc pr) (transpose xs)
                                                order = (map (\i-> (map (!!i) condE')) mapping, replicate (length mapping) priorE')

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
                                                 (condE'',priorE'') = stochmapOrder (condE',priorE') mapping priors
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

getScales :: Node -> Double -> (Double,Double)
getScales tree target = (x,x**n) where
        n = fromIntegral $ length $ leaves tree
        x = n `nthRootApprox` target 
