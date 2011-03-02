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

numberify :: ListAlignment -> [[Int]]
numberify aln = map nfy myseqs where
        myseqs = sequences aln
        nfy = numberMap 0
        numberMap i [] = []
        numberMap i ('-':xs) = -1 : numberMap i xs
        numberMap i (x:xs) = i : numberMap (i+1) xs

homDist :: ListAlignment -> ListAlignment -> Double
homDist aln1 aln2 = answer (diff3 pairs1 pairs2) where
        
        answer :: (Int,Int) -> Double
        answer (numPairs,numDiffs) = (fromIntegral numDiffs) / (fromIntegral numPairs)
        
        pairs :: [[Int]] -> [[[(Int,Int)]]]
        pairs [] = []
        pairs (x:[]) = []
        pairs (x:xs) = pairsXY (zip (repeat x) xs) : pairs (xs)

        pairsXY []  = []
        pairsXY ((a,b):xs) = zip a b : pairsXY xs

        pairs1 = pairs (numberify aln1)
        pairs2 = pairs (numberify aln2)

        addT (a,b) (c,d) = (a+c,b+d)

        diff2 :: [(Int,Int)] -> [(Int,Int)] -> (Int,Int)
        diff2 ((x1,x2):xs) y | x1 < 0 = addT (diff2 xs y) (0,0)
        diff2 x ((y1,y2):ys) | y1 < 0 = addT (diff2 x ys) (0,0)
        diff2 ((x1,x2):xs) ((y1,y2):ys) | x1 == y1 && x2==y2 = addT (diff2 xs ys) (1,0)
        diff2 ((x1,x2):xs) ((y1,y2):ys) | x1 == y1 = addT (diff2 xs ys) (1,1)
    --    diff2 (x:xs) (y:ys) = diff2 xs ys
        diff2 [] [] = (0,0)

        diff :: [[(Int,Int)]] -> [[(Int,Int)]] -> (Int,Int)
        diff [] [] = (0,0)
        diff (x:xs) (y:ys) = addT (diff2 x y) (diff xs ys)


        diff3 :: [[[(Int,Int)]]] -> [[[(Int,Int)]]] -> (Int,Int)
        diff3 [] [] = (0,0)
        diff3 (x:xs) (y:ys) = addT (diff x y) (diff3 xs ys)
