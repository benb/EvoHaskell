module Tree where
import Text.ParserCombinators.Parsec hiding (many, optional, (<|>))
import Control.Applicative hiding ((<|>))
import Text.ParserCombinators.Parsec.Expr
import Text.ParserCombinators.Parsec.Prim
import Numeric (readFloat, readHex, readSigned)
import Text.ParserCombinators.Parsec.Token
import Text.ParserCombinators.Parsec.Language
import qualified Data.HashMap as HM
import Data.List

data Node = Leaf {name :: String,distance :: Double } | INode Node Node Double | Tree Node Node deriving (Show,Eq) 

data PNode = PLeaf {pname :: String, pdistance :: Double } | PINode [PNode] Double | PTree [PNode] deriving Show


readNewickTree :: String -> Either ParseError PNode
readNewickTree = parse parseTree "" 

readBiNewickTree :: String -> Either ParseError Node
readBiNewickTree x = case (readNewickTree x) of  
                        Left err -> Left err
                        Right t -> Right (enforceBi t)

leaves :: Node -> [Node]
leaves = traverse []

traverse :: [Node] -> Node -> [Node]
traverse init (INode l r dist) = traverse (traverse init r) l 
traverse init (Tree l r) = traverse (traverse init r) l 
traverse init x = x : init

names :: Node -> [String]
names = map getName . leaves where
                getName (Leaf name dist) = name

enforceBi :: PNode -> Node
enforceBi (PLeaf name dist) = Leaf name dist
enforceBi (PINode (a:b:[]) dist) = INode (enforceBi a) (enforceBi b) dist
enforceBi (PINode (a:xs) dist) = INode (enforceBi a) (enforceBi $ PINode xs 0.0) dist 
enforceBi (PTree (a:b:[])) =  Tree (enforceBi a) (enforceBi b)
enforceBi (PTree (a:xs)) = Tree (enforceBi a) (enforceBi $ PINode xs 0.0)

-- http://evolution.genetics.washington.edu/phylip/newicktree.html :
-- "A name can be any string of printable characters except blanks, colons, semcolons, parentheses, and square brackets." 

nodeName :: Parser String
nodeName =  many1 (noneOf " \t\n\r:;()[]")

--from Real World Haskell
number :: Parser Double
number = do s <- getInput
            case readSigned readFloat s of
              [(n, s')] -> n <$ setInput s'
              _         -> fail "empty"
            

parseLeaf :: Parser PNode
parseLeaf = do name <- nodeName
               string ":"
               len <- number
               return (PLeaf name len)

parseINode :: Parser PNode
parseINode = do nodes <- parseGenINode 
                string ":"
                len <- number
                return (PINode nodes len)

parseTree :: Parser PNode
parseTree = do nodes <- parseGenINode
               option "" $ string ":0.0" --RAxML sticks branch lengths on the root...!
               string ";"
               return (PTree nodes)

parseGenINode :: Parser [PNode]
parseGenINode = do string "("
                   fst <- try (parseINode) <|> parseLeaf
                   string ","
                   remainder <- parseNodeList
                   string ")"
                   return (fst:remainder)

parseNodeList :: Parser [PNode]
parseNodeList = sepBy1 (try (parseINode) <|> parseLeaf) (char ',')

type Split = ([String],[String])

splits :: Node -> [Split]
splits = splits' [] []


splits' :: [Split] -> [String] -> Node -> [Split]
splits' init upper (Leaf name dist) =  init
splits' init upper (INode l r dist) = ((names l) ++ upper,(names r)): ((names r)++upper,(names l)) :  splits' (splits' init ((names r)++upper) l) ((names l)++upper) r
splits' init upper (Tree l r ) = ((names l),(names r)) :  splits' (splits' init ((names r)++upper) l) ((names l)++upper) r

splitsFor :: Node -> [String] -> HM.HashMap String Node
splitsFor node list = splitsFor' node HM.empty list

splitsFor' :: Node -> HM.HashMap String Node -> [String] -> HM.HashMap String Node
splitsFor' root startMap [] = startMap
splitsFor' root startMap remaining | (sort (names root)) == remaining = foldl' (\hash key -> HM.insert key root hash) startMap remaining

splitsFor' (Tree left right) startMap remaining = (go left right (go right left startMap)) where
                                                     go l r map = splitsFor' l map reduced where
                                                       reduced = remaining \\ (names r)
splitsFor' (INode left right dist) startMap remaining = (go left right (go right left startMap)) where
                                                     go l r map = splitsFor' l map reduced where
                                                       reduced = remaining \\ (names r)

splitsFor' leaf startMap remaining = error $ "Tree and Alignment are not congruent: Can't find leaf " ++ (head (remaining \\ (names leaf))) ++ " in tree"
