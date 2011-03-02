import System.Environment (getArgs)
import Alignment
import Text.Printf
main = do args <- getArgs
          (diff args) >>= putStrLn

diff (x:y:xs) = do a <- parseFastaFile x
                   b <- parseFastaFile y
                   return (printf "%.6f" (homDist a b))
diff x = return "Usage: homdist <fasta1> <fasta2>"
                        
       
