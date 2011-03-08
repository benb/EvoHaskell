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

gapPos :: Sequence -> [(Int,Int)]
gapPos s = gapPos' s [] Nothing 0

absoluteGapPos :: [(Int,Int)] -> [(Int,Int)]
absoluteGapPos s = map firstTwo (absoluteGapPos' s) where
                        firstTwo (a,b,c) = (a,b)

absoluteGapPos' :: [(Int,Int)] -> [(Int,Int,Int)]
absoluteGapPos' [] = [] 
absoluteGapPos' ((i,j):[]) = (i,j,0):[]
absoluteGapPos' ((i,j):xs) = (i+myoffset,j,myoffset) : agpTail where
                                agpTail = absoluteGapPos' xs
                                offset ((a,b,c):ys) = c + b
                                myoffset = offset agpTail


gapPos' :: Sequence -> [(Int,Int)] -> (Maybe Int) -> Int -> [(Int,Int)]
--end of sequence
gapPos' [] list Nothing pos = list 
gapPos' [] list (Just i) pos = (pos,i):list
--open a gap
gapPos' ('-':xs) list Nothing pos = gapPos' xs list (Just 1) (pos) 
--extend a gap
gapPos' ('-':xs) list (Just i) pos = gapPos' xs list (Just (i+1)) (pos)
--close a gap
gapPos' (x:xs) list (Just i) pos = gapPos' xs ((pos,i):list) Nothing (pos+1)
--no gap
gapPos' (x:xs) list Nothing pos = gapPos' xs list Nothing (pos+1)



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

