module Alignment where
import qualified Data.ByteString.Lazy.Char8 as L
import Control.Monad
import Data.List

appendString :: [(String,String)] -> String -> [(String, String)]
appendString old add = case  old of 
                (name,x):xs -> (name,(x++add)):xs

parseFasta :: [L.ByteString] -> [(String,String)] -> [(String,String)]
parseFasta [] old =  old
parseFasta bs old =  case L.unpack (L.take 1 (head bs)) of 
                      ['>'] -> parseFasta (tail bs) ((L.unpack (L.drop 1 (head bs)),"") : old)
                      _ -> parseFasta (tail bs) (appendString old (L.unpack (head bs))) 
 
parseFastaString :: L.ByteString -> ListAlignment

parseFastaString input = quickListAlignment names seqs where 
                                mydata = sortBy sortX (parseFasta (L.lines input) [])
                                sortX (a,b) (c,d) = compare a c
                                names = map fst mydata
                                seqs = map snd mydata 


parseFastaFile :: String -> IO ListAlignment
parseFastaFile name = parseFastaString `liftM` (L.readFile name)

type Name = String
type Sequence = [Char]
type Column = [Char]

data ListAlignment = ListAlignment {names ::  [Name],
                            sequences :: [Sequence],
                            columns :: [Column]} deriving Show

quickListAlignment :: [Name] -> [Sequence] -> ListAlignment
quickListAlignment names sequences = ListAlignment names sequences (transpose sequences)


toFasta :: ListAlignment -> [String]
toFasta aln = stringList where --foldl (++) "" stringList where 
                 stringList = map toSeqStr seqList
                 seqList = zip (names aln) (sequences aln)
                 toSeqStr :: (String,String) -> String
                 toSeqStr (name,seq) = ">" ++ name ++ "\n" ++ seq ++ "\n"

