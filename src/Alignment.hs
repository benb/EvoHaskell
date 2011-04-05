module Alignment where
import Tree hiding (names)
import qualified Tree
import qualified Data.ByteString.Lazy.Char8 as L
import Control.Monad
import Data.List
import Data.Char
import Debug.Trace
import qualified Data.HashMap as HM

appendString :: [(String,String)] -> String -> [(String, String)]
appendString old add = case  old of 
                (name,x):xs -> (name,x++add):xs

--Fasta format
parseFasta' :: [(String,String)] -> [L.ByteString] -> [(String,String)]
--this is an inefficient way to check that this alignment will be valid
parseFasta' ((name,""):(name2,seq2):(name3,seq3):xs) bs | (length seq2) /= (length seq3) = error $ "lengths of sequences " ++ name2 ++ " and " ++ name3 ++ " don't match"
parseFasta' old [] =  old
parseFasta' old bs =  case L.unpack (L.take 1 (head bs)) of 
                      ['>'] -> parseFasta' ((trim $ L.unpack (L.drop 1 (head bs)),"") : old) $ tail bs
                      _ -> parseFasta' (appendString old $ trim (L.unpack (head bs))) $ tail bs 

parseFasta :: Monad m => [L.ByteString] -> m [(String,String)]
parseFasta input = return $ parseFasta' [] input 

trim :: String -> String
trim = reverse . dropWhile isSpace . reverse . dropWhile isSpace 

 

parseAlignmentString :: Monad m =>  ([L.ByteString] -> m [(String,String)]) -> L.ByteString -> m ListAlignment

parseAlignmentString parser input = (liftM2 safeListAlignment names seqs) where 
                                mydata = (liftM $ sortBy sortX) (parser (L.lines input))
                                sortX (a,b) (c,d) = compare a c
                                names = (liftM $ map fst) mydata
                                seqs = (liftM $ map snd) mydata 


parseAlignmentFile :: Monad m => ([L.ByteString] -> m [(String,String)]) -> String -> IO (m ListAlignment)
parseAlignmentFile parser name = parseAlignmentString parser `liftM` (L.readFile name)

--"relaxed" phylip format

parsePhylip :: Monad m =>  [L.ByteString] -> m [(String,String)]
parsePhylip = parsePhylipHeader
parsePhylipHeader :: Monad m => [L.ByteString] -> m [(String,String)]
parsePhylipBody :: Int -> Int -> [(L.ByteString,L.ByteString)] -> [L.ByteString] -> [(String,String)]
parsePhylipBody' :: [L.ByteString] -> [(L.ByteString,L.ByteString)] -> [(L.ByteString,L.ByteString)] -> [(String,String)]


parsePhylipHeader (x:xs)  = (case (ans,nChar) of
                                (Just a,Just c) | find (\x-> length (snd x)/=c) a == Nothing -> return a
                                _ -> fail "Can't parse alignment") where
                                header = L.dropWhile (==' ') x
                                t1 = L.readInt header
                                nTaxa = fmap fst t1
                                t2 = fmap snd t1 >>= L.readInt . L.dropWhile  (==' ')
                                nChar = fmap fst t2
                                ans = fmap (\x -> parsePhylipBody x x [] xs) nTaxa

parsePhylipBody nTaxa 0 output (x:xs) = parsePhylipBody' xs [] $ reverse output
parsePhylipBody nTaxa remaining output (x:xs) = parsePhylipBody nTaxa (remaining-1) ((name,seq):output) xs where
                                                (name,remainder) = L.break (==' ') x
                                                seq = L.filter (/=' ') remainder
                                             
parsePhylipBody' [] top [] = map (\(x,y) -> (L.unpack x,L.unpack y)) $ reverse top 
parsePhylipBody' (x:xs) top [] = parsePhylipBody' (x:xs) [] $ reverse top
parsePhylipBody' (x:xs) top ((name,seq):ys) = parsePhylipBody' xs ((name,seq `L.append` (L.filter (/=' ') x)):top) ys

dropGaps :: ListAlignment -> [(String,String)]
dropGaps a = zip (names a) (map dropGap $ sequences a)

dropGap :: String -> String
dropGap xs = filter (not . isGapChar) xs

isGapChar :: Char -> Bool
isGapChar '-' = True
isGapChar '.' = True
isGapChar '~' = True
isGapChar x = False


type Name = String
type Sequence = [Char]
type Column = [Char]

data NumberedColumn = NumberedColumn {coldata::[(Char,Int)]} deriving (Eq)
instance Ord NumberedColumn where 
                       compare (NumberedColumn x) (NumberedColumn y) = compare' x y Nothing where
                                --maintain sequence ordering
                                compare' ((x,i):xs) ((y,j):ys) ans | (not $ isGapChar x) && (not $ isGapChar y) = compare i j

                                -- two all gap columns!?
                                compare' [] [] Nothing = EQ 

                                -- all pairs are gap-base or base-gap --> arbitrary
                                compare' [] [] (Just ans) = ans 

                                --first gap on left side -> set arbitrary answer and keep looking for base-base pairs
                                compare' ((gap,i):xs) ((y,j):ys) Nothing | (isGapChar gap) && (not $ isGapChar y) = compare' xs ys (Just GT)
                                compare' ((x,i):xs) ((gap,j):ys) Nothing | (isGapChar gap) && (not $ isGapChar x) = compare' xs ys (Just LT)
                                compare' ((gap1,i):xs) ((gap2,j):ys) ans | (isGapChar gap1) && (isGapChar gap2) =  compare' xs ys ans

                                --gap-base or base gap with existing ordering
                                compare' ((gap,i):xs) ((y,j):ys) ans | (isGapChar gap) = compare' xs ys ans
                                compare' ((x,i):xs) ((gap,j):ys) ans | (isGapChar gap) = compare' xs ys ans

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
gapPos' (gap:xs) list Nothing pos | isGapChar gap = gapPos' xs list (Just 1) (pos) 
--extend a gap
gapPos' (gap:xs) list (Just i) pos | isGapChar gap = gapPos' xs list (Just (i+1)) (pos)
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

safeListAlignment names sequences = removeAllGaps $ quickListAlignment names sequences

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
notAllGap (gap:[]) | isGapChar gap = False
notAllGap (gap:xs) | isGapChar gap = notAllGap xs
notAllGap (x:xs) = True

hasGap :: Column -> Bool
hasGap [] = False
hasGap (gap:xs) | isGapChar gap = True
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
        numberMap i (gap:xs) |isGapChar gap = -1 : numberMap i xs
        numberMap i (x:xs) = i : numberMap (i+1) xs

numberifyGap :: ListAlignment -> [[Int]]
numberifyGap aln = map nfy myseqs where
        myseqs = sequences aln
        nfy = numberMap 0
        numberMap i [] = []
        numberMap i (gap:xs) | isGapChar gap = (-(i+1)) : numberMap i xs
        numberMap i (x:xs) = i : numberMap (i+1) xs

numberifyGapTree :: Node -> ListAlignment -> [[(Int,Maybe Node)]]
numberifyGapTree tree aln = transpose $ nfy (columns aln) where
        nfy :: [Column]  -> [[(Int, Maybe Node)]]
        nfy colList = numberMap (map (\x->0) (head (columns aln))) colList
        numberMap :: [Int] -> [Column] -> [[(Int, Maybe Node)]]
        numberMap y [] = []
        numberMap y (x:xs) = (snd ans) : (numberMap (fst ans) xs) where
                              ans = numberMap' y x $ names aln
                              gapNums = unrootedSplitsFor tree gapNames
                              gapNames = map (\x-> fst x) $ filter (\t -> isGapChar (snd t)) $ zip (names aln) x
                              numberMap':: [Int] -> Column -> [String] -> ([Int],[(Int,Maybe Node)]) 
                              numberMap' [] [] [] = ([],[])

                              numberMap' (a:as) (gap:bs) (name:cs) | isGapChar gap = (a:(fst ans2),((-a-1,Just $ getNode name)):(snd ans2)) where
                                                                                             ans2 = numberMap' as bs cs
                              numberMap' (a:as) (b:bs) (name:cs) = (a+1:(fst ans2),(a,Nothing):(snd ans2)) where
                                                                                             ans2 = numberMap' as bs cs
                              getNode::String -> Node
                              getNode name = case (HM.lookup name gapNums) of 
                                Nothing -> error $ "Can't find gap for " ++ name
                                Just a -> a

compatible :: Node -> ListAlignment -> Bool
compatible tree aln = (sort $ Tree.names tree) == names aln


