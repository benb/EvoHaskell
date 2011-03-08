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
            

parseLeaf :: Parser Node
parseLeaf = do name <- nodeName
               string ":"
               len <- number
               return (Leaf name len)

parseINode :: Parser Node
parseINode = do (fst,snd) <- parseGenINode 
                string ":"
                len <- number
                return (INode fst snd len)

parseTree :: Parser Node
parseTree = do (fst,snd) <- parseGenINode
               string ";"
               return (Tree fst snd)

parseGenINode :: Parser (Node,Node)
parseGenINode = do string "("
                   fst <- try (parseINode) <|> parseLeaf
                   string ","
                   snd <- try (parseINode) <|> parseLeaf
                   string ")"
                   return  (fst,snd)

