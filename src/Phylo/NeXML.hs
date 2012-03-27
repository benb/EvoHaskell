module Phylo.NeXML where
import Phylo.Likelihood 
import Data.Maybe
import Control.Monad.State
import Control.Monad
import Debug.Trace


        
data FlatTree = FlatNode (Maybe String) Int | FlatEdge Int Int [Double]

nodes leaf@(DLeaf _ _ _ _ _ _) = [leaf]
nodes node@(DINode l r _ _ _) = node:((nodes l) ++ (nodes r))
nodes node@(DTree l m r _ _ _ _ ) = node:((nodes l) ++ (nodes m) ++ (nodes r))

setParent l (i,j,k) = (l,j,k)
xmlTree' :: DNode -> State (DNode,[String],Int) [String]


xmlTree' node@(DINode l r bl _ _ )  | trace (show node) True = do 
                                                               modify (setParent node) 
                                                               xmlTree' l
                                                               bll<-get
                                                               let bl = makeBranch bll l
                                                               modify (setParent node) 
                                                               xmlTree' r
                                                               brr<-get
                                                               let br = makeBranch brr r
                                                               (_,curr,id) <- get
                                                               let extra = map (\x -> x $ id+1) [br,bl]
                                                               let fin = (makeNode (id+1) Nothing):(extra ++ curr)
                                                               put (node, fin, id+1)
                                                               return fin



xmlTree' node@(DTree l m r _ _ _ _ )  = do 
                                       modify (setParent node) 
                                       xmlTree' l
                                       bll <- get
                                       let bl = makeBranch bll l
                                       modify (setParent node) 
                                       xmlTree' m
                                       bmm <- get
                                       let bm = makeBranch bmm m
                                       modify (setParent node) 
                                       xmlTree' r
                                       brr <- get
                                       let br = makeBranch brr r
                                       (_,curr,id) <- get
                                       let extra = map (\x -> x $ id+1) [br,bm,bl]
                                       let fin = (makeNode (id+1) Nothing):(extra ++ curr)
                                       put (node,fin,id+1)
                                       return $ fin


xmlTree' node@(DLeaf name bl _ _ _ _ )  = do 
                                         (parent,curr,id) <- get
                                         put (node,(makeNode (id+1) (Just name)):curr,id+1)
                                         return (name:curr)

makeBranch (_,_,to) (DINode _ _ bl _ _ ) from = makeBranch' from to bl
makeBranch (_,_,to) (DLeaf _ bl _ _ _ _ ) from = makeBranch' from to bl
makeBranch' from to bl = "<edge source=\"" ++ (show from) ++ "\" target=\"" ++ (show to) ++ "\" length=\""  ++ (show $ head bl) ++"\" />"

makeNode :: Int -> Maybe String -> String
makeNode id (Just name) = "<node id=\"" ++ (show id) ++ "\" label=" ++ (show name) ++ "/>"
makeNode id Nothing = "<node id=\"" ++ (show id) ++ "\" />"

xmlTree :: DNode -> String
xmlTree node = unlines $ xmlHeader ++ (evalState (xmlTree' node) (node,[],0)) ++ xmlTail

xmlHeader = ["<?xml version=\"1.0\" encoding=\"UTF-8\"?>","<nex:nexml xmlns:nex=\"http://www.nexml.org/2009\" xmlns=\"http://www.nexml.org/2009\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:xsd=\"http://www.w3.org/2001/XMLSchema#\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" generator=\"org.nexml.model.impl.DocumentImpl\" version=\"0.9\">","  <otus id=\"otus1\"/>","  <trees id=\"trees2\" otus=\"otus1\">","    <tree about=\"#tree3\" id=\"tree3\" label=\"Tree1\" xsi:type=\"nex:FloatTree\">"]
xmlTail = ["    </tree>", "  </trees>", "</nex:nexml>"]
