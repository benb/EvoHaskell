module Phylo.Alignment.Parsers where
import qualified Data.ByteString.Lazy.Char8 as L
import Data.List
import Data.Char

--This parser has been briefly tested to work with the "relaxed" phylip format
--It is intended to handle interleaved and sequential versions, but relies on spaces between 
--The sequence name and the sequence data
--It also won't work if sequence names are repeated
parsePhylip :: Monad m =>  [L.ByteString] -> m [(String,String)]
parsePhylip = parsePhylipHeader
parsePhylipHeader :: Monad m => [L.ByteString] -> m [(String,String)]
parsePhylipBody :: Bool -> Int -> Int -> [(L.ByteString,L.ByteString)] -> [L.ByteString] -> [(String,String)]
parsePhylipBody' :: [L.ByteString] -> [(L.ByteString,L.ByteString)] -> [(L.ByteString,L.ByteString)] -> [(String,String)]


parsePhylipHeader (x:xs)  = (case (ans,nChar) of
                                (Just a,Just c) | find (\x-> length (snd x)/=c) a == Nothing -> return a
                                _ -> fail "Can't parse alignment") where
                                header = L.dropWhile (==' ') x
                                t1 = L.readInt header
                                nTaxa = fmap fst t1
                                t2 = fmap snd t1 >>= L.readInt . L.dropWhile  (==' ')
                                nChar = fmap fst t2
                                filteredxs = filter (\x->L.empty /= (L.filter (/=' ') x)) xs
                                ans = fmap (\x -> parsePhylipBody False x x [] filteredxs) nTaxa

--parsePhylipBody nTaxa remaining output xs | trace ((show nTaxa) ++ " " ++ (show remaining) ++ " " ++ (show $  length xs) ++ (show (take 1 xs)) ++ (show output)) False  = undefined
parsePhylipBody True nTaxa 0 output xs =  parsePhylipBody' [] [] $ reverse finoutput where
                                                                finoutput = (name,seq):(tail output)
                                                                (name,seq1) = head output
                                                                seq=appendall seq1 xs
                                                                appendall a [] = a
                                                                appendall a (y:ys) = appendall (a `L.append` (L.filter (/=' ') y)) ys
                                                                
parsePhylipBody False nTaxa 0 output xs =  parsePhylipBody' xs [] $ reverse output
--new sequence
parsePhylipBody sequential nTaxa remaining output (x:xs) | (L.head x) /= ' ' && notseen x output = parsePhylipBody sequential nTaxa (remaining-1) ((name,seq):output) xs where
                                                                                                                (name,remainder) = L.break (==' ') x
                                                                                                                seq = L.filter (/=' ') remainder

--seen before, append to previous seq
parsePhylipBody sequential nTaxa remaining output (x:xs) | (L.head x) /= ' ' = parsePhylipBody sequential nTaxa remaining output (droppedx:xs) where
                                                                               (_,droppedx)=L.break(==' ') x

--no name so seen before
parsePhylipBody sequential nTaxa remaining output (x:xs)         = parsePhylipBody True nTaxa remaining ((name,seq):(tail output)) xs where
                                                                                (name,seq1) = head output
                                                                                seq=seq1 `L.append` (L.filter(/=' ') x)
notseen:: L.ByteString -> [(L.ByteString,L.ByteString)] -> Bool
notseen line previous = not $ any (name==) names where
                                names = map fst previous
                                (name,_) = L.break (==' ') line
                                             
--params are (remaining lines) (output stack 1) (output stack 2)
--we are finished:
parsePhylipBody' [] top [] = map (\(x,y) -> (L.unpack x,L.unpack y)) $ reverse top 
parsePhylipBody' [] [] out = map (\(x,y) -> (L.unpack x,L.unpack y)) $ out

--loop next set of sequences
parsePhylipBody' (x:xs) top [] = parsePhylipBody' (x:xs) [] $ reverse top

--append x to sequence seq and push onto top
parsePhylipBody' (x:xs) top ((name,seq):ys) = parsePhylipBody' xs ((name,seq `L.append` (cleanup x)):top) ys where
                                                                            -- remove (repeated) name if present 
                                                                cleanup str | L.head str /=' ' = cleanup $ snd $ L.break (==' ') str
                                                                            | otherwise = L.filter(/=' ') str




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

appendString :: [(String,String)] -> String -> [(String, String)]
appendString old add = case  old of 
                (name,x):xs -> (name,x++add):xs

 
