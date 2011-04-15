import System.Environment (getArgs)
import System.Console.GetOpt
import Phylo.Alignment
import Data.List
import Control.Monad
import Phylo.Tree
import Phylo.Data
import Phylo.Likelihood


data Flag = AlignmentFile String | TreeFile String
options = [ Option ['a'] ["alignment"] (ReqArg AlignmentFile "FILE") "Alignment",
            Option ['t'] ["tree"] (ReqArg TreeFile "FILE") "Tree" ]

main = do args <- getArgs
          (aln,tree) <- case getOpt Permute options args of 
                         (((AlignmentFile aln):(TreeFile tre):[]),[],[]) -> do aln <- parseAlignmentFile parseFasta aln
                                                                               tree <- (liftM readBiNewickTree) (readFile tre)
                                                                               return (aln,tree)
                         (_,_,msgs) -> error $ concat msgs
          case (aln,tree) of 
                (Just a,Right t)-> putStrLn $ "WAG lkl:" ++ (show $ quickLkl a t wagPi wagS)
                (_,_) -> error "Can't parse something"

 

