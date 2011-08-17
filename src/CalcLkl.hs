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


data Flag = AlignmentFile String | TreeFile String | Alpha String | OptAlpha | OptThmm
options = [ Option ['a'] ["alignment"] (ReqArg AlignmentFile "FILE") "Alignment",
            Option ['t'] ["tree"] (ReqArg TreeFile "FILE") "Tree",
            Option ['g'] ["gamma"] (ReqArg Alpha "DOUBLE") "Use Gamma Model" ,
            Option ['o'] ["opt-gamma"] (NoArg OptAlpha) "Optimise Gamma Model",
            Option ['m'] ["opt-thmm"] (NoArg OptThmm) "Optimise THMM Model" ]

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
                                                        f = optGammaF 4 a t wagPi wagS
                                                        alpha = goldenSection 0.001 0.02 99.0 (\x -> -(f x))
                                                        lkl = f alpha
                (Just a,Right t,(OptThmm):[])-> putStrLn $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show sigma) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        t2 = structDataN 5 AminoAcid (pAlignment a) t                                                                                                                                                                       
                                                        piF = fromList $ scaledAAFrequencies a
                                                        model = thmmModel 5 piF wagS
                                                        (optTree,[priorZero,logalpha,logsigma],_) = optParamsAndBL model t2 [0.1,0.1,0.2] [1.0] [Just 0.01,Nothing,Nothing] [Just 0.99,Nothing,Nothing] 0.01
                                                        (alpha,sigma) = (exp logalpha,exp logsigma)
                                                        lkl = -(logLikelihood optTree)
                (_,_,_) -> error "Can't parse something"

 

