import System.Environment (getArgs)
import Alignment
import Alignment.Dist
import Text.Printf
import System.Console.GetOpt

main = do args <- getArgs
          let (dist,files) = if ((head args)=="-g") 
                then (homGapDist,(drop 1 args))
                else (homDist,args)
          (diff dist files) >>= putStrLn

diff dist (x:y:xs) = do a <- parseFastaFile x
                        b <- parseFastaFile y
                        return (printf "%.6f" (dist a b))
diff dist x = return "Usage: homdist <fasta1> <fasta2>"
                        
       
