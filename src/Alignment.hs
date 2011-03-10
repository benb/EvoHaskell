module Alignment where
import Tree hiding (names)
import qualified Tree
import qualified Data.ByteString.Lazy.Char8 as L
import Control.Monad
import Data.List
import qualified Data.HashMap as HM

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

data NumberedColumn = NumberedColumn {coldata::[(Char,Int)]} deriving (Eq)
instance Ord NumberedColumn where 
                       compare (NumberedColumn x) (NumberedColumn y) = compare' x y Nothing where
                                --maintain sequence ordering
                                compare' ((x,i):xs) ((y,j):ys) ans | x/='-' && y/='-' = compare i j

                                -- two all gap columns!?
                                compare' [] [] Nothing = EQ 

                                -- all pairs are gap-base or base-gap --> arbitrary
                                compare' [] [] (Just ans) = ans 

                                --first gap on left side -> set arbitrary answer and keep looking for base-base pairs
                                compare' (('-',i):xs) ((y,j):ys) Nothing | y/='-' = compare' xs ys (Just GT)
                                compare' ((x,i):xs) (('-',j):ys) Nothing | x/='-'= compare' xs ys (Just LT)
                                compare' (('-',i):xs) (('-',j):ys) ans = compare' xs ys ans

                                --gap-base or base gap with existing ordering
                                compare' (('-',i):xs) ((y,j):ys) ans = compare' xs ys ans
                                compare' ((x,i):xs) (('-',j):ys) ans = compare' xs ys ans

sortAlignment :: ListAlignment -> ListAlignment
sortAlignment (ListAlignment names seqs cols) = ListAlignment names (transpose ans) ans where
                                                  numbers = transpose $ numberifyBasic $ ListAlignment names seqs cols
                                                  numbCols = map NumberedColumn $ map (\(a,b)-> zip a b) $ zip cols numbers
                                                  reordered = sort numbCols
                                                  ans = map (map fst) (map coldata reordered)



                                                



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
instance Eq ListAlignment where 
        (==) a b = (names a) == (names b) && (sequences a) == (sequences b) --assume cols are ok


quickListAlignment :: [Name] -> [Sequence] -> ListAlignment
quickListAlignment names sequences = ListAlignment names sequences (transpose sequences)

fromColumnListAlignment :: [Name] -> [Column] -> ListAlignment
fromColumnListAlignment names cols = ListAlignment names (transpose cols) cols


toFasta :: ListAlignment -> [String]
toFasta aln = stringList where --foldl (++) "" stringList where 
                 stringList = map toSeqStr seqList
                 seqList = zip (names aln) (sequences aln)
                 toSeqStr :: (String,String) -> String
                 toSeqStr (name,seq) = ">" ++ name ++ "\n" ++ seq ++ "\n"

removeAllGaps :: ListAlignment -> ListAlignment
removeAllGaps (ListAlignment names seqs cols) = fromColumnListAlignment names (removeAllGaps' cols)

removeAllGaps' :: [Column] -> [Column]
removeAllGaps' = filter notAllGap where 
 
notAllGap :: Column -> Bool
notAllGap ('-':[]) = False
notAllGap ('-':xs) = notAllGap xs
notAllGap (x:xs) = True

hasGap :: Column -> Bool
hasGap [] = False
hasGap ('-':xs) = True
hasGap (x:xs) = hasGap xs


--orderGaps :: ListAlignment -> ListAlignment
--orderGaps (ListAlignment names seqs cols) = fromColumnListAlignment names (orderGaps' cols)


--orderGaps' :: [Column] -> [Column] -> [Column]
--orderGaps' [] x = x:[]
--orderGaps' (x:xs) [] | hasGap x = orderGaps' xs x:[]
--orderGaps' (x:xs) [] = orderGaps' xs []
--orderGaps (x:xs) (y:ys) | canPushTogether x y 
                        
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

numberifyGapTree :: Node -> ListAlignment -> [[(Int,Maybe Node)]]
numberifyGapTree tree aln = transpose $ nfy (columns aln) where
        nfy :: [Column]  -> [[(Int, Maybe Node)]]
        nfy colList = numberMap (map (\x->0) (head (columns aln))) colList
        numberMap :: [Int] -> [Column] -> [[(Int, Maybe Node)]]
        numberMap y [] = []
        numberMap y (x:xs) = (snd ans) : (numberMap (fst ans) xs) where
                              ans = numberMap' y x $ names aln
                              gapNums = splitsFor tree gapNames
                              gapNames = map (\x-> fst x) $ filter (\t -> (snd t)=='-') $ zip (names aln) x
                              numberMap':: [Int] -> Column -> [String] -> ([Int],[(Int,Maybe Node)]) 
                              numberMap' [] [] [] = ([],[])

                              numberMap' (a:as) ('-':bs) (name:cs) = (a:(fst ans2),((-a-1,Just $ getNode name)):(snd ans2)) where
                                                                                             ans2 = numberMap' as bs cs
                              numberMap' (a:as) (b:bs) (name:cs) = (a+1:(fst ans2),(a,Nothing):(snd ans2)) where
                                                                                             ans2 = numberMap' as bs cs
                              getNode::String -> Node
                              getNode name = case (HM.lookup name gapNums) of 
                                Nothing -> error $ "Can't find gap for " ++ name
                                Just a -> a

compatible :: Node -> ListAlignment -> Bool
compatible tree aln = (sort $ Tree.names tree) == (names aln)


