import System.Environment (getArgs)
import Alignment
main = do args <- getArgs
          (diff args) >>= print 

diff (x:y:xs) = do a <- parseFastaFile x
                   b <- parseFastaFile y
                   return (homDist a b)
       
