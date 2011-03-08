module Tree where
import Control.Applicative
import Text.ParserCombinators.Parsec hiding (many, optional, (<|>))
import Text.ParserCombinators.Parsec.Expr
import Numeric (readFloat, readHex, readSigned)
import Text.ParserCombinators.Parsec.Token
import Text.ParserCombinators.Parsec.Language

data Node = Leaf {name :: String,distance :: Double } | INode Node Node Double | Tree Node Node deriving Show 

data PNode = PLeaf {pname :: String, pdistance :: Double } | PINode [PNode] Double | PTree [PNode] deriving Show

leaves :: Node -> [Node]
leaves (Leaf name dist)  = [Leaf name dist]
leaves (INode fst snd dist) = (leaves fst) ++ (leaves snd)
leaves (Tree fst snd) = (leaves fst) ++ (leaves snd)

enforceBi :: PNode -> Node
enforceBi (PLeaf name dist) = Leaf name dist
enforceBi (PINode (a:b:[]) dist) = INode (enforceBi a) (enforceBi b) dist
enforceBi (PINode (a:xs) dist) = INode (enforceBi a) (enforceBi $ PINode xs 0.0) dist 
enforceBi (PTree (a:b:[])) =  Tree (enforceBi a) (enforceBi b)
enforceBi (PTree (a:xs)) = Tree (enforceBi a) (enforceBi $ PINode xs 0.0)

--parseTree :: Parser Tree
--parseTree = do {
--                char '('
--                left <- parseNode
--                char ','
--                right <- parseNode
--                char ')'
--                char ';'
--                return Tree left right
--} 
--
-- http://evolution.genetics.washington.edu/phylip/newicktree.html :
-- "A name can be any string of printable characters except blanks, colons, semcolons, parentheses, and square brackets." 

nodeName :: Parser String
nodeName =  many1 (noneOf " \t\n\r:;()[]")

--from Real World Haskell
number :: Parser Double
number = do s <- getInput
            case readSigned readFloat s of
              [(n, s')] -> n <$ setInput s'
              _         -> empty
            

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


