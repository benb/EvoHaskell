import System.Environment (getArgs)
import Alignment
main = do
         (parseFastaFile "in.fa") >>= mapM putStr . toFasta 
       
