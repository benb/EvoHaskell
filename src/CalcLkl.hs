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


data Flag = AlignmentFile String | TreeFile String | Alpha String | OptAlpha | OptThmm | OptThmmP | OptThmm2 String
options = [ Option ['a'] ["alignment"] (ReqArg AlignmentFile "FILE") "Alignment",
            Option ['t'] ["tree"] (ReqArg TreeFile "FILE") "Tree",
            Option ['g'] ["gamma"] (ReqArg Alpha "DOUBLE") "Use Gamma Model" ,
            Option ['o'] ["opt-gamma"] (NoArg OptAlpha) "Optimise Gamma Model",
            Option ['m'] ["opt-thmm"] (NoArg OptThmm) "Optimise THMM Model" ,
            Option ['n'] ["opt-thmmplus"] (NoArg OptThmmP) "Optimise THMM Model",
            Option ['p'] ["opt-thmm2"] (ReqArg OptThmm2 "splitsStr") "2 state THMM"]

main = do args <- getArgs
          (aln,tree,remainOpts) <- case getOpt Permute options args of 
                         (((AlignmentFile aln):(TreeFile tre):xs),[],[]) -> do aln <- parseAlignmentFile parseUniversal aln
                                                                               tree <- (liftM readBiNewickTree) (readFile tre)
                                                                               return (aln,tree,xs)
                         (_,_,msgs) -> error $ concat msgs
          case (aln,tree,remainOpts) of 
                (Just a,Right t,[])-> putStrLn $ "WAG lkl:" ++ (show $ quickLkl a t wagPi wagS)
                (Just a,Right t,(Alpha val):[])-> putStrLn $ "WAG Gamma lkl:" ++ (show $ quickGamma 4 (read val) a t wagPi wagS)
                (Just a,Right t,(OptAlpha):[])-> putStrLn $ "Opt Alpha: " ++ (show alpha) ++ " " ++ (show lkl) where
                                                        t2 = structDataN 1 AminoAcid (pAlignment a) t
                                                        piF = fromList $ scaledAAFrequencies a
                                                        model = gammaModel 4 wagS piF
                                                        (optTree,[logalpha],_) = optParamsAndBL model t2 [0.5] [0.25,0.25,0.25,0.25] [Just 0.001] [Nothing] 0.01
                                                        alpha = exp logalpha
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(OptThmm):[])-> putStrLn $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show sigma) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        t2 = structDataN 5 AminoAcid (pAlignment a) t                                                                                                                                                                       
                                                        piF = fromList $ scaledAAFrequencies a
                                                        model = thmmModel 5 wagS piF
                                                        (optTree,[priorZero,alpha,sigma],_) = optParamsAndBL model t2 [0.1,0.5,1.0] [1.0] [Just 0.01,Just 0.001, Just 0.00] [Just 0.99,Nothing,Nothing] 0.01
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(OptThmmP):[])-> putStrLn $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show optTree) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        piF = fromList $ scaledAAFrequencies a
                                                        t2 = addModelFx (structDataN 5 AminoAcid (pAlignment a) t) (gammaModel 4 wagS piF [0.5]) $ flatPriors 4                                                                                                                                                                
                                                        startBL = getBL t2
                                                        newBL = map (\x -> x:0.0:[]) startBL
                                                        t3 = fst $ setBLX' 0 newBL t2
                                                        model = thmmPerBranchModel 5 wagS piF
                                                        (optTree,[priorZero,alpha],_) = optParamsAndBL model t3 [0.1,0.5] [1.0] [Just 0.01,Just 0.001] [Just 0.99,Nothing] 0.01
                                                        bls = getBL t3
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(OptThmm2 spl):[])-> putStrLn $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show optTree) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        piF = fromList $ scaledAAFrequencies a
                                                        t2 = addModelFx (structDataN 5 AminoAcid (pAlignment a) t) (gammaModel 4 wagS piF [0.5]) $ flatPriors 4                                                                                                                                                                
                                                        startBL = getBL t2
                                                        newBL = map (\x -> x:0.0:[]) startBL
                                                        t3 = fst $ setBLX' 0 newBL t2
                                                        goodNodes = map (filter (/=',')) $ groupBy (\x y -> y /= ',') spl
                                                        mapped = makeMapping (\(x,y) -> if (((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))) then (trace "OK0" 0) else (trace "OK1" 1)) t3
                                                        model | traceShow t3 True = thmmPerBranchModel 5 wagS piF 
                                                        (optTree,[priorZero,alpha,sigma0,sigma1]) = optBSParams model t3 [0.1,0.5,1.0,2.0] 2 mapped [1.0] [Just 0.01,Just 0.001,Just 0.0,Just 0.0] [Just 0.99,Nothing,Nothing,Nothing] 0.01
                                                       -- optTree = optBLD0 $ (addModelFx t3 (model [priorZero,alpha])) bls
                                                        bls = getBL t3
                                                        lkl = logLikelihood optTree
                (_,_,_) -> error "Can't parse something"

 

