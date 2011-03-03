module Alignment.Dist where
import Alignment

numberifyBasic :: ListAlignment -> [[Int]]
numberifyBasic aln = map nfy myseqs where
        myseqs = sequences aln
        nfy = numberMap 0
        numberMap i [] = []
        numberMap i ('-':xs) = -1 : numberMap i xs
        numberMap i (x:xs) = i : numberMap (i+1) xs

numberifyGap :: ListAlignment -> [[Int]]
numberifyGap aln = map nfy myseqs where
        myseqs = sequences aln
        nfy = numberMap 0
        numberMap i [] = []
        numberMap i ('-':xs) = (-(i+1)) : numberMap i xs
        numberMap i (x:xs) = i : numberMap (i+1) xs

homDist = genDist numberifyBasic
homGapDist = genDist numberifyGap

genDist :: (ListAlignment -> [[Int]]) -> ListAlignment -> ListAlignment -> Double
genDist numF aln1 aln2 = answer (diff3 pairs1 pairs2) where
        
        pairs1 = pairs (numF aln1)
        pairs2 = pairs (numF aln2)

        answer :: (Int,Int) -> Double
        answer (numPairs,numDiffs) = (fromIntegral numDiffs) / (fromIntegral numPairs)
        
        pairs :: [[Int]] -> [[[(Int,Int)]]]
        pairs [] = []
        pairs (x:[]) = []
        pairs (x:xs) = pairsXY (zip (repeat x) xs) : pairs (xs)

        pairsXY []  = []
        pairsXY ((a,b):xs) = zip a b : pairsXY xs

        addT (a,b) (c,d) = (a+c,b+d)

        diff2 :: [(Int,Int)] -> [(Int,Int)] -> (Int,Int)
        diff2 ((x1,x2):xs) y | x1 < 0 = diff2 xs y
        diff2 x ((y1,y2):ys) | y1 < 0 = diff2 x ys
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
