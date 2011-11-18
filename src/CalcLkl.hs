{-# LANGUAGE ForeignFunctionInterface #-}
import Foreign
import Foreign.C.Types
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
foreign import ccall "stochmap.h CalculateAndWrite"
     c_CalculateAndWrite :: CInt -> CInt -> CInt -> CInt -> CInt -> Ptr (Ptr (Ptr (Ptr CInt))) -> Ptr (Ptr CInt) -> Ptr CInt -> Ptr CInt -> Ptr (Ptr (Ptr (Ptr (Ptr CDouble)))) -> Ptr (Ptr (Ptr CDouble)) -> Ptr (Ptr CDouble) -> Ptr (Ptr CDouble) -> Ptr CDouble -> Ptr CDouble -> Ptr CFile -> IO ()



data Flag = NumCats String | AlignmentFile String | TreeFile String | Alpha String | OptAlpha | OptThmm | OptThmmP | OptThmm2 String | OptSim1 | OptSim2 String | OptSim0 String | Seed String
options = [ Option ['a'] ["alignment"] (ReqArg AlignmentFile "FILE") "Alignment",
            Option ['t'] ["tree"] (ReqArg TreeFile "FILE") "Tree",
            Option ['s'] ["seed"] (ReqArg Seed "INTEGER") "Seed",
            Option [] ["num-cats"] (ReqArg NumCats "NUMGAMMACATS") "Number of Gamma Rate Categories",
            Option ['g'] ["gamma"] (ReqArg Alpha "DOUBLE") "Use Gamma Model" ,
            Option ['o'] ["opt-gamma"] (NoArg OptAlpha) "Optimise Gamma Model",
            Option ['m'] ["opt-thmm"] (NoArg OptThmm) "Optimise THMM Model" ,
            Option ['n'] ["opt-thmmplus"] (NoArg OptThmmP) "Optimise THMM Model",
            Option ['p'] ["opt-thmm2"] (ReqArg OptThmm2 "splitsStr") "2 state THMM",
            Option [] ["sim0"] (ReqArg OptSim0 "PARAMS") "Simulate WAG",
            Option [] ["sim1"] (NoArg OptSim1) "Simulation 1",
            Option [] ["sim2"] (ReqArg OptSim2 "PARAMS") "Simulation 2"]

main = do args <- getArgs
          (aln,tree,stdGen,remainOpts) <- case getOpt Permute options args of 
                         (((AlignmentFile aln):(TreeFile tre):(Seed seed):xs),[],[]) -> do aln <- parseAlignmentFile parseUniversal aln
                                                                                           tree <- (liftM readBiNewickTree) (readFile tre)
                                                                                           let stdGen = mkStdGen $ read seed
                                                                                           return (aln,tree,stdGen,xs)
                         (((AlignmentFile aln):(TreeFile tre):xs),[],[]) -> do aln <- parseAlignmentFile parseUniversal aln
                                                                               tree <- (liftM readBiNewickTree) (readFile tre)
                                                                               stdGen <- getStdGen
                                                                               return (aln,tree,stdGen,xs)
                         (_,_,msgs) -> error $ concat msgs
          let (cats,remainOpts') = case remainOpts of 
                (NumCats n):rem -> (read n,rem)
                rem -> (4,rem)
          case (aln,tree,remainOpts') of 
                (Just a,Right t,[])-> putStrLn $ "WAG lkl:" ++ (show $ quickLkl a t wagPi wagS)
                (Just a,Right t,(Alpha val):[])-> putStrLn $ "WAG Gamma lkl:" ++ (show $ quickGamma cats (read val) a t wagPi wagS)
                (Just a,Right t,(OptAlpha):[])-> putStrLn $ "Opt Alpha: " ++ (show alpha) ++ " " ++ (show lkl) where
                                                        t2 = structDataN 1 AminoAcid (pAlignment a) t
                                                        piF = fromList $ scaledAAFrequencies a
                                                        model = gammaModel cats wagS piF
                                                        (optTree,[logalpha],_) = optParamsAndBL model t2 [0.5] [0.25,0.25,0.25,0.25] [Just 0.001] [Nothing] 0.01
                                                        alpha = exp logalpha
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(OptThmm):[])-> putStrLn $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show sigma) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        t2 = structDataN 5 AminoAcid (pAlignment a) t                                                                                                                                                                       
                                                        piF = fromList $ scaledAAFrequencies a
                                                        model = thmmModel 5 wagS piF
                                                        (optTree,[priorZero,alpha,sigma],_) = optParamsAndBL model t2 [0.1,0.5,1.0] [1.0] [Just 0.01,Just 0.001, Just 0.00] [Just 0.99,Nothing,Nothing] 0.01
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(ThmmStochmap params ):[])-> let alpha:sigma:priorZero:[] = map read $ words params
                                                             c_CalculateAndWrite nSite nState nBranch nProc nCols scales lMat multiplicites sitemap partials qset sitelikes pi_i branchLengths mixProbs where
                                                               pAln = pAlignment a
                                                               (nSiteRaw,nColsRaw,multiplicitesRaw) = getAlnData pAln where
                                                                  getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts)
                                                               model = thmmModel 5 wagS piF
                                                               t2 = structDataN 5 AminoAcid (pAln) t                                                                                                                                                                       
                                                               piF = fromList $ scaledAAFrequencies a
                                                               nSite <- newArray nSiteRaw
                                                               nCols <- newArray nColsRaw
                                                               multiplicites <- newArray multiplicitesRaw
                                                               lMat <- newArray2 $ makeIntraLMat 5 20
                                                               nProc = length model
                                                               pi_i <- newArray piF
                                                               branchLengths <- newArray $ getBL t2
                                                               mixProbs <- newArray $ getPriors t2
                                                        (optTree,[priorZero,alpha,sigma],_) = optParamsAndBL model t2 [0.1,0.5,1.0] [1.0] [Just 0.01,Just 0.001, Just 0.00] [Just 0.99,Nothing,Nothing] 0.01
                (Just a,Right t,(OptThmmP):[])-> putStrLn $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show optTree) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        piF = fromList $ scaledAAFrequencies a
                                                        t2 = addModelFx (structDataN 5 AminoAcid (pAlignment a) t) (gammaModel cats wagS piF [0.5]) $ flatPriors cats                                                                                                                                                                
                                                        startBL = getBL t2
                                                        newBL = map (\x -> x:0.0:[]) startBL
                                                        t3 = fst $ setBLX' 0 newBL t2
                                                        model = thmmPerBranchModel 5 wagS piF
                                                        (optTree,[priorZero,alpha],_) = optParamsAndBL model t3 [0.1,0.5] [1.0] [Just 0.01,Just 0.001] [Just 0.99,Nothing] 0.01
                                                        bls = getBL t3
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(OptThmm2 spl):[])-> putStrLn $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show optTree) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        piF = fromList $ scaledAAFrequencies a
                                                        t2 = addModelFx (structDataN 5 AminoAcid (pAlignment a) t) (gammaModel cats wagS piF [0.5]) $ flatPriors cats                                                                                                                                                             
                                                        startBL = getBL t2
                                                        newBL = map (\x -> x:0.0:[]) startBL
                                                        t3 = fst $ setBLX' 0 newBL t2
                                                        goodNodes = map (filter (/=',')) $ groupBy (\x y -> y /= ',') spl
                                                        mapped = makeMapping (\(x,y) -> if (((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))) then (trace "OK0" 0) else (trace "OK1" 1)) t3
                                                        model | traceShow t3 True = thmmPerBranchModel 5 wagS piF 
                                                        (optTree,[priorZero,alpha,sigma0,sigma1]) = optBSParams model t3 [0.1,0.5,2.0,2.0] 2 mapped [1.0] [Just 0.01,Just 0.001,Just 0.0,Just 0.0] [Just 0.99,Nothing,Nothing,Nothing] 0.0001
                                                       -- optTree = optBLD0 $ (addModelFx t3 (model [priorZero,alpha])) bls
                                                        bls = getBL t3
                                                        lkl = logLikelihood optTree
                (Just a,Right t,OptSim1:[]) -> putStrLn $ concat $ toFasta $ simulateSequences AminoAcid 5 stdGen 349 t2 where
                                                        alpha = 0.8335109218715334
                                                        priorZero = 0.1907077998691572
                                                        sigma0 = 15.631076379213784
                                                        sigma1 = 1.0836675920110694e-8
                                                        piF = fromList $ scaledAAFrequencies a
                                                        t2 = addModelFx (setBLMapped 1 (dummyTree (structDataN 5 AminoAcid (pAlignment a) t)) mapped ) (thmmPerBranchModel 5 cpRevS cpRevPi [priorZero,alpha]) [1.0]
                                                        goodNodes = ["E_Nosloc","E_Enccun","E_Gluple","A_Aerper","A_Metbar"] 
                                                        mapped = makeMapping (\(x,y) -> if (((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))) then ([sigma0]) else ([sigma1])) t2
                (Just a,Right t,(OptSim0 params):[]) -> putStrLn $ concat $ toFasta $ simulateSequences AminoAcid 5 stdGen 349 t2 where
                                                        alpha:sigma0:sigma1:priorZero:[] = map (read) $ words params
                                                        piF = fromList $ scaledAAFrequencies a
                                                        t2 = addModelFx (setBLMapped 1 (dummyTree (structDataN 5 AminoAcid (pAlignment a) t)) mapped ) (thmmPerBranchModel 5 wagS piF [priorZero,alpha]) [1.0]
                                                        goodNodes = ["E_Nosloc","E_Enccun","E_Gluple","A_Aerper","A_Metbar"] 
                                                        mapped = makeMapping (\(x,y) -> if (((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))) then ([sigma0]) else ([sigma1])) t2
                (Just a,Right t,OptSim2 params:[]) -> do (lgInnerS,lgInnerPi) <- parsePamlDatIO "lgInner"
                                                         (lgOuterS,lgOuterPi) <- parsePamlDatIO "lgOuter"
                                                         putStrLn $ compute (lgInnerS,lgInnerPi) (lgOuterS,lgOuterPi) where
                                                             compute (s0,pi0) (s1,pi1) = concat $ toFasta $ simulateSequences AminoAcid 5 stdGen 349 t3 where
                                                                alpha:sigma0:sigma1:priorZero:[] = map (read) $ words params
                                                                piF = fromList $ scaledAAFrequencies a
                                                                (model0,pi_0) = thmmPerBranchModel 5 s0 pi0 [priorZero,alpha]
                                                                (model1,pi_1) = thmmPerBranchModel 5 s1 pi1 [priorZero,alpha]
                                                                t2 = addModelFx (setBLMapped 1 (dummyTree (structDataN 5 AminoAcid (pAlignment a) t)) mapped ) (thmmPerBranchModel 5 cpRevS cpRevPi [priorZero,alpha]) [1.0]
                                                                tfFunc x y = ((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))
                                                                goodNodes = ["E_Nosloc","E_Enccun","E_Gluple","A_Aerper","A_Metbar"] 
                                                                mapped = makeMapping (\(x,y) -> if (tfFunc x y)  then ([sigma0]) else ([sigma1])) t2
                                                                mappedModels = makeMapping (\(x,y) -> if (tfFunc x y) then model1 else model0) t2
                                                                t3 = restructDataMapped t2 mappedModels [1.0] pi_0
                (_,_,_) -> error "Can't parse something"

newArray2 :: [[Storable]] -> IO Ptr Ptr Storable
newArray2 array = do aList <- mapM newArray array
                     array <- newArray aList
                     return array

makeIntraLMat -> Int -> Int -> [[Int]]
makeIntraLMat nClass alphabet = makeIntraLMat' nClass alphabet (nClass*alphabet - 1) [] 
makeIntraLMat' nClass alphabet -1 xs = xs
makeIntraLMat' nClass alphabet pos xs = makeIntraLMat' nClass alphabet (pos-1) ((getRow pos):xs)  where
        getRow i = (replicate (chunk * 20) 0) ++ allOne ++ (replicate (nClass - chunk -1) * 20) where
                chunk = i `div` (nClass * alphabet)
                allOne = replicate 20 1

