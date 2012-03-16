module Phylo.Likelihood where
import Phylo.Alignment
import Phylo.Tree
import Phylo.Matrix
import Phylo.Opt
import Phylo.NLOpt
import Numeric.LinearAlgebra ((<>),(<.>),scale,mul,constant,diag)
import Data.Packed.Matrix
import Data.Packed.Vector
import Data.List
import Data.Maybe
import Statistics.Math
import Numeric.GSL.Distribution.Continuous
import Numeric.GSL.Minimization
import Control.Exception as E
import Control.Parallel
import Control.Parallel.Strategies
import Debug.Trace
import Phylo.Data
import System.Random
import Numeric.GSL.Differentiation


-- |Combine two sets of partial likelihoods along two branches
combinePartial :: Matrix Double -> Matrix Double -> Matrix Double
combinePartial = mul 

-- |Combine two lists of partial likelihoods (for mixture models)
allCombine :: [Matrix Double] -> [Matrix Double] -> [Matrix Double]
allCombine (left:lefts) (right:rights) = (combinePartial left right):(allCombine lefts rights)
allCombine [] [] = []

rawPartialLikelihoodCalc :: Matrix Double -> Matrix Double -> Matrix Double
rawPartialLikelihoodCalc pT start = pT <> start

allPLC :: [Matrix Double] -> [Matrix Double] -> [Matrix Double]
allPLC (pT:pTs) (start:starts) = (rawPartialLikelihoodCalc pT start):(allPLC pTs starts)
allPLC [] [] = []


sumLikelihoods :: [Int] -> [Double] -> Double
sumLikelihoods counts likelihoods = foldr foldF 0 $ zip (map fromIntegral counts) likelihoods where
                                        foldF :: (Double,Double) -> Double -> Double 
                                        foldF (i,y) x = x + (i * (log y))
                                        

logLikelihoodModel :: DNode -> Double
logLikelihoodModel root@(DTree _ _ _ _ pC _ _) = sumLikelihoods pC $ likelihoods root where

concatNewLine :: [[String]] -> String
concatNewLine = intercalate "\n" . map (intercalate "\n")


posteriorTipsCSV :: Int -> PatternAlignment -> DNode -> String
posteriorTipsCSV i aln j = intercalate "\n" $ (posteriorTipsCSVX aln) $ posteriorTips i aln j

posteriorTipsCSVX aln = map (posteriorTipsCSVX' aln)

posteriorTipsCSVX' :: PatternAlignment ->  [(Int,String,DNode,[Double])] -> String
posteriorTipsCSVX' (PatternAlignment names seqs _ _ _)  list = intercalate "\n" $ map (\(x,x',y,z) -> name ++ "," ++  (show x) ++ "," ++ (show x') ++"," ++ [z]  ++"," ++ (intercalate "," (map show y))) $ zip4 [0..] realseqpos (transpose post) alignedseq where
                         post = map (\(x,y,z,lkl)->lkl) list
                         name = head $ map (\(x,y,z,lkl) -> y) list
                         id = fromJust $ findIndex (==name) names
                         alignedseq = seqs !! id
                         realseqpos = reverse $ foldl incPos [] alignedseq
                         incPos [] '-' = [-1]
                         incPos [] base = [0]
                         incPos (x:xs) '-' = (x:x:xs) 
                         incPos (x:xs) base = (x+1:x:xs)


posteriorTipsCSV' :: (Int,String,DNode,[Double]) -> String
posteriorTipsCSV' (i,j,k,lkl) = j ++ "," ++ (show i) ++ "," ++ (intercalate "," $ map show lkl)

posteriorTips :: Int -> PatternAlignment -> DNode -> [[(Int,String,DNode,[Double])]]
posteriorTips numModels aln root@(DTree l m r pLs pC priors pis) = map normaliseL2 $ posteriorTipsLkl numModels aln root where
                                                                --   getLkl (i,j,k,lkl) = lkl
                                                                --   getLklList = (map . map) getLkl
normaliseL2 :: [(Int,String,DNode,[Double])] -> [(Int,String,DNode,[Double])]
normaliseL2 list = replaced where
                    unzipped = unzip4 list
                    (_,_,_,lkls) = unzipped
                    probs = normaliseL2' lkls
                    replaced = map (\((i,j,k,lkl),probs) -> (i,j,k,probs)) $ zip list probs

normaliseL2' :: [[Double]] -> [[Double]]
normaliseL2' list = transpose $ normaliseL2'' $ transpose list
normaliseL2'' = map normaliseL
normaliseL :: [Double] -> [Double]                                                                        
normaliseL xs = normalise' xs (sum xs)
normalise' xs tot = map (/tot) xs

posteriorTipsLkl :: Int -> PatternAlignment -> DNode -> [[(Int,String,DNode,[Double])]]
posteriorTipsLkl numModels aln root@(DTree _ _ _ _ _ _ _) = map (map (\(i,j,k)->(i,j,k,likelihoodsPerSite aln k)))$ posteriorTipsList numModels root

posteriorTipsList :: Int -> DNode -> [[(Int,String,DNode)]]
posteriorTipsList numModels (DLeaf name dist seq partial model pl) = [map getLeaf [0..(numModels-1)]] where
                                                                        getLeaf i = (i,name,(DLeaf name dist seq partial' model (calcLeafPL partial' dist model))) where
                                                                                partial' = blankMat partial numModels i

posteriorTipsList numModels (DINode l r bl model pl) = (map (map getINodeR) $ posteriorTipsList numModels r)++ (map (map getINodeL) $ posteriorTipsList numModels l) where
                                                           getINodeL (i,name,node) = (i,name,DINode node r bl model (calcPL node r bl model))
                                                           getINodeR (i,name,node) = (i,name,DINode l node bl model (calcPL l node bl model))

posteriorTipsList numModels (DTree l m r pLs model priors pis) = (map (map getINodeR) $ posteriorTipsList numModels r) ++ (map (map getINodeM) $ posteriorTipsList numModels m) ++ (map (map getINodeL) $ posteriorTipsList numModels l) where
                                                           getINodeL (i,name,node) = (i,name,DTree node m r (calcRootPL node m r) model priors pis)
                                                           getINodeM (i,name,node) = (i,name,DTree l node r (calcRootPL l node r) model priors pis)
                                                           getINodeR (i,name,node) = (i,name,DTree l m node (calcRootPL l m node) model priors pis)
                                                                
                                                           


blankMat :: Matrix Double -> Int -> Int -> Matrix Double
blankMat matrix numClass except = fromRows filtRows where
                                        blocksize = (rows matrix) `div` numClass
                                        lower = except * blocksize
                                        upper = lower + blocksize
                                        myRows = toRows matrix
                                        zeroVect = mapVector (\x->0.0) (head myRows)
                                        filtRows = map (\(i,j) -> if (j >= lower && j < upper) then i else zeroVect) (zip myRows [0..])
                                        

likelihoods :: DNode -> [Double]
likelihoods (DTree _ _ _ pLs _ priors pis) = map summarise (transpose likelihoodss) where
                                                summarise lkl = sum $ zipWith (*) lkl priors
                                                likelihoodss = map toList $ map (\(pi,pL) -> pi <> pL) $ zip pis pLs

likelihoodsPerSite aln tree = likelihoodss' where
                                likelihoodss = likelihoods tree
                                mapping = mapBack aln
                                likelihoodss' = map (likelihoodss !!) mapping
                                                                

logLikelihood = logLikelihoodModel

data PatternAlignment = PatternAlignment {names :: [String], seqs::[String], columns::[String], patterns::[String], counts::[Int]}

pAlignment (ListAlignment names seqs columns) = PatternAlignment names seqs columns patterns counts where
                                                   (patterns,counts) = unzip $ map (\x -> (head x,length x)) $ group $ sort columns
mapBack :: PatternAlignment -> [Int]
mapBack (PatternAlignment _ _ columns patterns _) = map (\x->fromJust $ findIndex (==x) patterns) columns

calcPL :: DNode -> DNode -> [Double] -> [BranchModel] -> [Matrix Double]
calcPL left right dist model = allPLC (map (\x-> fst $ x dist) model) $ allCombine (getPL left) (getPL right) 
calcLeafPL :: Matrix Double -> [Double] -> [BranchModel] -> [Matrix Double]
calcLeafPL tips dist model = map (\pT -> rawPartialLikelihoodCalc pT tips) $ (map (\x -> fst $ x dist) model)
calcRootPL left middle right = allCombine (getPL middle) $ allCombine (getPL left) (getPL right)

{-- gives [[[Matrix Double]]] where [branch [proc [left/right Matrix]]] --}
{-- need to map from [branch [l/r [proc Matrix]]] --}
toPBEList ((f,s,_,_,_):xs) = (map (\(i,j) -> [i,j]) $ zip f s):(toPBEList xs)
toPBEList [] = []

toPBELengths ((_,_,dist,_,_):xs) = dist:(toPBELengths xs)
toPBELengths [] = []

toPBESplits = map (\(_,_,_,a,_) -> a)
toPBEQ :: [([Matrix Double],[Matrix Double],Double,([String],[String]),[Matrix Double])] -> [[Matrix Double]]
toPBEQ = map (\(_,_,_,_,a) -> a)

--getQ node = concat $ getQ' node [] 
getQ = getQ'
getQ' (DLeaf _ bl  _  _ model _)   = (map (\x-> snd $ x bl) model) 
getQ' (DINode l r bl model _) = map (\x-> snd $ x bl) model


getPartialBranchEnds :: DNode -> [([Matrix Double],[Matrix Double],Double,([String],[String]),[Matrix Double])]
getPartialBranchEnds tree = getPartialBranchEnds' $ allRootings tree
getPartialBranchEnds' (x:xs) = (getPartialBranchEndsLeaves'' x) ++ (getPartialBranchEndsX' xs) -- first tree is special case, only do leaves
getPartialBranchEndsX' (x:xs) = (getPartialBranchEnds'' x) ++ (getPartialBranchEndsX' xs)
getPartialBranchEndsX' [] = []
getPartialBranchEnds'' (DTree left middle right pLS _ _ _) = (branchPL right,noright,head $ getBL right,splits,getQ right): remainder where
                                                                noright = allCombine (getPL left) (getPL middle) 
                                                                remainder = (calc left middle) ++ (calc middle left)
                                                                calc node other = case (node,other) of 
                                                                                       (l@(DLeaf _ bl _ leafPL _ _),o) -> [(branchPL l,allCombine (getPL o) (getPL right),head $ getBL l,((descendents o) ++ (descendents right), descendents l),getQ l)]
                                                                                       (_,o) -> []
                                                                splits = ((descendents left) ++ (descendents middle), (descendents right))

getPartialBranchEndsLeaves'' (DTree left middle right pLS _ _ _) = remainder where
                                                                noright = allCombine (getPL left) (getPL middle) 
                                                                remainder = (calc right middle left) ++ (calc left middle right) ++ (calc middle left right)
                                                                calc node other other'= case (node,other,other') of 
                                                                                       (l@(DLeaf _ bl _ leafPL _ _),o,o') -> [(branchPL l,allCombine (getPL o) (getPL o'),head $ getBL l,((descendents o) ++ (descendents o'), descendents l),getQ l)]
                                                                                       (_,_,_) -> []


branchPL :: DNode -> [Matrix Double]
branchPL (DLeaf _ _ _ leafPL _ endPL) = map (\x->leafPL) endPL
branchPL (DINode l r _ _ endPL) = map (\(x,y)-> combinePartial x y) (zip (getPL l) (getPL r))

restructData node model priors pi = restructDataMapped node (\x -> model) priors pi

restructDataMapped :: DNode -> (DNode -> [BranchModel]) -> [Double] -> [Vector Double] -> DNode 
restructDataMapped leaf@(DLeaf name dist sequence partial _ _) model priors pi  = DLeaf name dist sequence partial myModel partial' where
                                                                               partial' = calcLeafPL partial dist myModel
                                                                               myModel = model leaf

restructDataMapped inode@(DINode left right dist _ _ ) model priors pi= newleft `par` newright `par` DINode newleft newright dist myModel partial where
                                                            partial = calcPL newleft newright dist myModel
                                                            myModel = model inode
                                                            newleft = restructDataMapped left model priors pi
                                                            newright = restructDataMapped right model priors pi
restructDataMapped tree@(DTree left middle right _ pc _ _ ) model priors pi = newleft `par` newright `par` newmiddle `par` DTree newleft newmiddle newright partial pc priors pi where 
                                            partial = calcRootPL newleft newmiddle newright
                                            newleft = restructDataMapped left model priors pi
                                            newright = restructDataMapped right model priors pi
                                            newmiddle = restructDataMapped middle model priors pi  



structDataN :: Int -> SeqDataType -> PatternAlignment -> Node -> NNode

structDataN hiddenClasses seqDataType pAln node = structDataN' hiddenClasses seqDataType pAln (transpose $ patterns pAln) node 
structDataN' hiddenClasses seqDataType (PatternAlignment names seqs columns patterns counts) transPat (Leaf name dist) = ans where
                                                                                                                          partial = getPartial hiddenClasses seqDataType sequence 
                                                                                                                          sequence = snd $ fromJust $ find (\(n,_) -> n==name) $ zip names transPat 
                                                                                                                          ans = NLeaf name [dist] sequence partial 


structDataN' hiddenClasses seqDataType pA transPat (INode c1 c2 dist) = left `par` right `par` NINode left right [dist] where
                                                                              left = structDataN' hiddenClasses seqDataType pA transPat c1
                                                                              right = structDataN' hiddenClasses seqDataType pA transPat c2

structDataN' hiddenClasses seqDataType pA transPat (Tree (INode l r dist) c2 ) = left `par` middle `par` right `par` NTree left middle right $ counts pA where
                                                                                  left = structDataN' hiddenClasses seqDataType pA transPat l
                                                                                  right = structDataN' hiddenClasses seqDataType pA transPat r
                                                                                  middle = structDataN' hiddenClasses seqDataType pA transPat $ addDist dist c2
                                                                                  addDist d (Leaf name dist) = Leaf name (dist+d)
                                                                                  addDist d (INode l r dist) = INode l r (dist+d)

structDataN' hiddenClasses seqDataType pA transPat (Tree c2 (INode l r dist)) = structDataN' hiddenClasses seqDataType pA transPat (Tree (INode l r dist) c2)

structDataN' hiddenClasses seqDataType pA transPat (Tree (Leaf name dist) (Leaf name2 dist2))  = error "Can't handle tree with only 2 leaves"


structData hiddenClasses seqDataType pAln model transPat node = addModel (structDataN hiddenClasses seqDataType pAln node) model 



addModelNNode :: NNode -> [BranchModel] -> [Double] -> [Vector Double] -> DNode
addModelNNode (NLeaf name dist sequence partial) model _ _  = DLeaf name dist sequence partial model partial' where
                                                                                              partial' = calcLeafPL partial dist model

addModelNNode (NINode c1 c2 dist) model priors pi =  left `par` right `par` DINode left right dist model myPL where
                                                  left = addModelNNode c1 model priors pi
                                                  right = addModelNNode c2 model priors pi
                                                  myPL = calcPL left right dist model 
                                                                            
addModelNNode (NTree c1 c2 c3 pat) model priors pi = left `par` right `par` middle `par` DTree left middle right myPL pat priors pi where
                                  left = addModelNNode c1 model priors pi
                                  middle = addModelNNode c2 model priors pi
                                  right= addModelNNode c3 model priors pi
                                  myPL = calcRootPL left middle right

data SeqDataType = AminoAcid | Nucleotide
getPartial :: Int -> SeqDataType -> String -> Matrix Double
getPartial a b c = trans $ fromRows $ getPartial' a b c
getPartial' :: Int -> SeqDataType -> String -> [Vector Double]
getPartial' _ _ [] = []
getPartial' classes AminoAcid (x:xs) = (aaPartial classes x):(getPartial' classes AminoAcid xs)
getPartial' classes Nucleotide (x:xs) = error "Nucleotides unimplemented"

aaOrder = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
aaPartial classes x | isGapChar x = buildVector (20*classes) (\i->1.0) 
                    | otherwise = case findIndex (==x) aaOrder of 
                                       Just index -> buildVector (20*classes) (\i -> if i`mod`20==index then 1.0 else 0.0) where
                                       Nothing -> error $ "character " ++ [x] ++ "is not an amino acid or gap"

                                 


quickGamma' numCat tree priors pi s aln alpha = logLikelihood dataTree  where
                                                dataTree = structData 1 AminoAcid pAln branchModel transpats tree priors (repeat pi)
                                                branchModel = (map (\r->(\x->(scaledpT r eigenS x,mat))) rates) 
                                                rates = gamma numCat alpha
                                                transpats = transpose $ patterns pAln
                                                (eigenS,mat) = quickEigen' pi s
                                                pAln = pAlignment aln
                                                patcounts = counts pAln
                                                
                                                        
quickEigen pi s = eigQ (normQ (makeQ s pi) pi) pi
quickEigen' pi s = (eigQ mat pi,mat) where
                        mat = normQ (makeQ s pi) pi

standardpT = scaledpT 1.0
scaledpT scale eigenS params = pT eigenS $ (head params) * scale


optGammaF :: Int -> ListAlignment -> Node -> Vector Double -> Matrix Double -> (Double -> Double)
optGammaF numCat aln tree pi s = quickGamma' numCat tree priors pi s aln where
                                                priors = take numCat $ repeat (1.0/(fromIntegral numCat))


--flatFullPi numCat pi = Data.Packed.Vector.mapVector (/(fromIntegral numCat)) $ Data.Packed.Vector.join $replicate numCat pi
flatFullPi numCat pi = makeFullPi (replicate numCat (1.0/(fromIntegral numCat))) pi
makeFullPi priors pi = Data.Packed.Vector.join $ map (\p->Data.Packed.Vector.mapVector (*p) pi) priors
flatPriors numCat = replicate numCat (1.0/ (fromIntegral numCat))

gammaMix numCat alpha (u,lambda,u') = map (\s -> (u,scale s lambda,u')) scales where
                           scales = gamma numCat alpha

gamma :: Int -> Double -> [Double]
gamma numCat shape | shape > 50000.0 = gamma numCat 50000.0 --work around gsl convergence errors (who wants >100 in the gamma dist anyway?)
                   | otherwise       = map rK' [0..(numCat-1)] where
                        alpha = shape
                        beta = shape
                        factor = alpha/beta*(fromIntegral numCat)

                        freqK = map freqKf [0..(numCat-1)]
                        freqKf i = incompleteGamma (alpha + 1.0) (beta * (gammaInvCDF (((fromIntegral i)+1.0)/(fromIntegral numCat))))

                        rK' 0 = (head freqK) * factor
                        rK' n | n < numCat = ((freqK!!n) - (freqK!!(n-1))) * factor
                              | otherwise = factor * (1.0-(freqK!!(numCat-2)))
                        gammaInvCDF p = (density_1p ChiSq UppInv (2.0*alpha) (1.0-p)) / (2.0*beta) -- (UppInv mysteriously has more stable convergence)


type BranchModel = [Double] -> (Matrix Double,Matrix Double) -- first matrix exponentiated for BL, second not

type OptModel = DNode -> [Double] -> Double
data DNode = DLeaf {dName :: String,dDistance :: [Double],sequence::String,tipLikelihoods::Matrix Double,p::[BranchModel],partialLikelihoods::[Matrix Double]} |
        DINode DNode DNode [Double] [BranchModel] [Matrix Double] | 
        DTree DNode DNode DNode [Matrix Double] [Int] [Double] [(Vector Double)]

data NNode = NLeaf String [Double] String (Matrix Double) | NINode NNode NNode [Double] | NTree NNode NNode NNode [Int]
data DataModel = DataModel {dTree::DNode, patterncounts::[Int], priors::[Double], pis::[Vector Double]}

getLeaves node = getLeaves' node []
getLeaves' :: DNode -> [DNode] -> [DNode]
getLeaves' (DTree l m r _ _ _ _) xs = getLeaves' r $ getLeaves' m $ getLeaves' l xs
getLeaves' (DINode l r _ _ _) xs = getLeaves' l $ getLeaves' r xs
getLeaves' leaf xs = (leaf:xs)

instance Show DNode where
        show (DTree l m r _ _ _ _) = "("++(show l)++","++(show m)++","++(show r)++");"
        show (DINode l r [d] _ _) = "(" ++ (show l)++","++(show r)++"):"++(show d)
        show (DINode l r d _ _) = "(" ++ (show l)++","++(show r)++"):"++(show d)
        show (DLeaf name [d] _ _ _ _) = name++":"++(show d)
        show (DLeaf name d _ _ _ _) = name++":"++(show d)


class AddTree a where
        addModel :: (AddTree a) => a -> [BranchModel] -> [Double] -> [Vector Double] -> DNode
        addModelFx :: (AddTree a) => a -> ([BranchModel],[Vector Double]) -> [Double] -> DNode
        addModelFx t (bm,pi) prior = addModel t bm prior pi
instance AddTree NNode where
        addModel w x y z = ans2 where
                       ans = addModelNNode w x y z
                       ans2 = seq (getPL ans) ans
instance AddTree DNode where
        addModel w x y z = ans2 where 
                           ans = restructData w x y z
                           ans2 = seq (getPL ans) ans



getPL (DLeaf _ _ _ _ _ lkl) = lkl
getPL (DINode _ _ _ _ lkl) = lkl
getPL (DTree _ _ _ lkl _ _ _) = lkl

getDist (DLeaf _ dist _ _ _ _) = dist
getDist (DINode _ _ dist _ _) = dist
getDist root = []

getModel (DLeaf _ _ _ _ model _) = model
getModel (DINode _ _ _ model _) = model
getModel root = error "No model for root"


goLeft (DTree (DINode l r dist model pL) middle right _ pC priors pi)  = DTree l r r' rootpL pC priors pi where
        rootpL = calcRootPL l r r'
        r' = DINode middle right dist model pL' 
        pL' = calcPL middle right dist model

goRight (DTree left middle (DINode l r dist model pL) _ pC priors pi) = DTree l' l r rootpL pC priors pi where
        rootpL = calcRootPL l' l r
        l' = DINode left middle dist model pL'
        pL' = calcPL left middle dist model

canGoLeft (DTree (DINode _ _ _ _ _) _ _ _ _ _ _ ) = True
canGoLeft _ = False
canGoMid (DTree _ (DINode _ _ _ _ _) _ _ _ _ _ ) = True
canGoMid _ = False

optBLD = optBLDx optLeftBL
optBLD0 = optBLDx (optLeftBLx 0)

optBLDFull0 :: DNode -> DNode
optBLDFull0 = optBLDFull  (optLeftBLx 0) 
optBLDFull f tree = optBLDFull' (logLikelihood tree) 1E-05 f tree
optBLDFull' old cutoff f tree = case (new - old) of 
        x | x > cutoff -> optBLDFull' new cutoff f newtree
        otherwise      -> newtree 
        where newtree = optBLDx f tree
              new = logLikelihood newtree
        
                                

optBLDx f tree = t3 where
        t1 = swapLeftMiddle $ f $ if (canGoLeft tree) then goRight $ optBL' f (goLeft tree) else tree
        t2 = swapLeftRight $ f $ swapLeftMiddle $ f $ if (canGoLeft t1) then goRight $ optBL' f (goLeft t1) else t1
        t3 = swapLeftRight $ f $ if (canGoLeft t2) then goRight $ optBL' f (goLeft t2) else t2

optBL' f tree = ans where
        leftOptd = if (canGoLeft tree) then goRight $ optBL' f (goLeft tree) else tree
        ansL = swapLeftMiddle $ f leftOptd
        midOptd = if (canGoLeft ansL) then goRight $ optBL' f (goLeft ansL) else ansL
        ans = swapLeftMiddle $ f midOptd

allRootings tree = ans where
        midTree = swapLeftMiddle tree
        rightTree = swapLeftRight tree
        leftRootings = if (canGoLeft tree) then (tree' : allLM (tree')) else [] 
        midRootings = if (canGoLeft midTree) then (midTree' : allLM (midTree')) else [] 
        rightRootings = if (canGoLeft rightTree) then (rightTree' : allLM (rightTree')) else [] 
        tree' = goLeft tree
        midTree' = goLeft midTree
        rightTree' = goLeft rightTree
        ans  = tree : (leftRootings ++ midRootings ++ rightRootings)


allLM tree = a ++ b where
               a = allLeftRootings tree
               b = allMidRootings tree

allMidRootings tree | canGoMid tree = t2 : (allLM t2) where
                                      t2 = goLeft $ swapLeftMiddle tree
allMidRootings tree = []

allLeftRootings tree | canGoLeft tree = t2 : (allLM t2) where
                                        t2 = goLeft tree
allLeftRootings tree = []

descendents (DTree l m r _ _ _ _) = (descendents l) ++ (descendents m ) ++ (descendents r)
descendents (DLeaf name _ _ _ _ _) = [name] 
descendents (DINode l r _ _ _) = (descendents l) ++ (descendents r)

nonRootNodes (DTree l m r _ _ _ _) = (nonRootNodes l) ++ (nonRootNodes m) ++ (nonRootNodes r)
nonRootNodes node@(DINode l r _ _ _) = node : ((nonRootNodes l) ++ (nonRootNodes r))
nonRootNodes node@(DLeaf _ _ _ _ _ _) = [node]

--setLeftBL (DTree l m r _ _ _ _ ) x | trace ((show x ) ++ " (" ++ (descendents l) ++ " : " ++ (descendents m) ++ " : " ++ (descendents r) ++ ")") False = undefined

setLeftBL (DTree (DLeaf name dist seq tip model _) m r _ pC priors pi) x = (DTree newleft m r pl' pC priors pi) where
        partial' = calcLeafPL tip x model
        newleft = DLeaf name x seq tip model partial'
        pl' = calcRootPL newleft m r 


setLeftBL (DTree (DINode l1 r1 dist model _) m r _ pC priors pi) x = (DTree newleft m r pl' pC priors pi) where
        partial' = calcPL l1 r1 x model
        newleft = DINode l1 r1 x model partial'
        pl' = calcRootPL newleft m r 

optLeftBL tree  = traceShow ("OptBL " ++ (show best) ++ "->" ++ (show $ logLikelihood $ setLeftBL tree best)) $ setLeftBL tree best where
                        (DTree l _ _ _ _ _ _ ) = tree
                        best = case l of 
                                (DINode _ _ [param] _ _) -> (safeGoldenSection 0.01 1E-6 10.0 (invert $ (\x -> logLikelihood $ setLeftBL tree [x]))) : []
                                (DLeaf _ [param] _ _ _ _) -> (safeGoldenSection 0.01 1E-6 10.0 (invert $ (\x -> logLikelihood $ setLeftBL tree [x]))) : []
                                (DLeaf _ params _ _ _ _) -> fst $ maximize NMSimplex2 1E-4 1000 (map (\i->0.05) params) (boundedBLfunc (\x -> logLikelihood $ setLeftBL tree x)) params    
                                (DINode _ _ params _ _) -> fst $ maximize NMSimplex2 1E-4 1000 (map (\i->0.05) params) (boundedBLfunc (\x -> logLikelihood $ setLeftBL tree x)) params    

getLeftBL (DTree (DINode _ _ x _ _) _ _ _ _ _ _) = x
getLeftBL (DTree (DLeaf _ x _ _ _ _) _ _ _ _ _ _) = x

optLeftBLx val tree  = traceShow ("OptBL " ++ (show startBL) ++ " -> " ++ (show best) ++ " -> " ++ (show $ logLikelihood $ setLeftBL tree best)) $ bestTree where
                                startBL = getLeftBL tree
                                startLL = logLikelihood tree
                                b = goldenSection 0.01 1E-6 10.0 (invert $ (\x -> logLikelihood $ setLeftBL tree (replace val [x] startBL)))
                                best = replace val [b] startBL
                                bestTree' = setLeftBL tree best 
                                --protect from golden section screw-ups
                                bestTree = case (logLikelihood bestTree') of 
                                                x | x < startLL -> tree
                                                  | otherwise -> bestTree'

boundedBLfunc = boundedFunction (-1E20) (repeat $ Just 0.0) (repeat Nothing)

swapLeftMiddle (DTree l m r pL pC priors pi) = DTree m l r pL pC priors pi
swapLeftRight (DTree l m r pL pC priors pi) = DTree r m l pL pC priors pi


getSplits node = map invertSplit $ getSplits' node [] where
                        alldescendents = sort $ descendents node
                        invertSplit s = ((sort s), filter (\x -> not $ x `elem` s) alldescendents) 

getSplits' (DTree l m r _ _ _ _ ) splits = (getSplits' l (getSplits' m (getSplits' r splits)))
getSplits' (DINode l r _ _ _ ) splits = (descendents l ++ descendents r):(getSplits' l (getSplits' r splits))
getSplits' (DLeaf name _ _ _ _ _) splits = [name]:splits

getBL node = concat $ getBL' node []
getPriors (DTree _ _ _ _ _ priors _ ) = priors

getBL' (DLeaf _ bl _ _ _ _) bls = bl : bls
getBL' (DINode l r bl _ _) bls = bl : (getBL' l (getBL' r bls))
getBL' (DTree l m r _ _ _ _ ) bls = (getBL' l (getBL' m (getBL' r bls)))

getLinearMap mapping (DTree l m r _ _ _ _ ) ints = (getLinearMap mapping l (getLinearMap mapping m (getLinearMap mapping r ints)))
getLinearMap mapping node@(DINode l r _ _ _) ints = (mapping node) : (getLinearMap mapping l (getLinearMap mapping r ints))
getLinearMap mapping leaf ints = (mapping leaf) : ints


{-setBL :: [Double] -> DNode -> DNode -}
{-setBL bls node = fst $ setBL' reformattedBLs node where-}
        {-oldBLs = getBL' node []-}
        {-reformattedBLs = reformat bls oldBLs-}
        {-reformat blList (x:xs) = (taken : reformat remainder xs) where-}
                {-(taken,remainder) = splitAt (length x) blList-}

setBL' x = setBLX' 0 x

setBLMapped i (DTree l m r _ pC priors pi) mapping = DTree left middle right partial pC priors pi where
        left = setBLMapped i l mapping
        middle = setBLMapped i m mapping
        right = setBLMapped i r mapping
        partial = calcRootPL left middle right

setBLMapped i (DINode l r blstart mats pl) mapping = DINode left right newbl mats partial where
        left = setBLMapped i l mapping
        right = setBLMapped i r mapping
        partial = calcPL left right (newbl) mats
        initBL = take i blstart
        endBL = mapping (DINode left right blstart mats pl)
        newbl = initBL ++ endBL

setBLMapped i (DLeaf a blstart b tips model pl) mapping = DLeaf a newbl b tips model partial where
        partial = calcLeafPL tips newbl model
        initBL = take i blstart
        endBL = mapping (DLeaf a blstart b tips model pl)
        newbl = initBL ++ endBL



setBLX' i bls (DTree l m r pl pC priors pi) = ((DTree left middle right partial pC priors pi), remainder3) where
                        (left,remainder) = setBLX' i bls l
                        (middle,remainder2) = setBLX' i remainder m
                        (right,remainder3) = setBLX' i remainder2 r
                        partial = calcRootPL left middle right

setBLX' i (bl2:bls) (DINode l r blstart mats pl) = ((DINode left right newbl mats partial), remainder2) where
                             (left,remainder) = setBLX' i bls l
                             (right,remainder2) = setBLX' i remainder r
                             partial = calcPL left right (newbl) mats
                             newbl = replace i bl2 blstart

setBLX' i (bl2:bls) (DLeaf a (blstart) b tips model _) = ((DLeaf a newbl b tips model pl),bls) where
                                      pl = calcLeafPL tips newbl model 
                                      newbl = replace i bl2 blstart

setBLX' _ [] _ = error "Insuffiencent branch lengths specified"


replace 0 [x] (s:ss) = (x:ss)
replace i x start =  (take (i) start) ++ x ++ (drop (i + (length x)) start)

listpartition (c:counts) start ans = listpartition counts b (a:ans) where
                                 (a,b) = splitAt c start
listpartition [] [] ans = reverse ans
listpartition [] start ans = listpartition [] [] (start:ans)


setBLflat i node = setBLflat' i node $ map (\x -> (length x) - i) (getBL' node [])
setBLflat' i node counts bls = setBLX' i (listpartition counts bls []) node
                        

zeroQMat size = diag $ constant 0.0 size

qdLkl numCat priors dataType aln tree modelF params = logLikelihood $ addModelFx (structDataN numCat dataType (pAlignment aln) tree) (modelF params) priors

type ModelF = [Double] -> ([BranchModel],[Vector Double])
basicModel :: Matrix Double -> Vector Double -> ModelF
basicModel s pi _ = ([model],[pi]) where
                    model = (\x -> (standardpT e x,m))
                    (e,m) = quickEigen' pi s

quickLkl aln tree pi s = qdLkl 1 [1.0] AminoAcid aln tree (basicModel s pi) []
quickGamma numCat alpha aln tree pi s = qdLkl 1 (flatPriors numCat) AminoAcid aln tree (gammaModel numCat s pi) [alpha]

gammaModel numCat s pi (alpha:xs) = (models,replicate numCat pi) where
                                models = map (\r -> (\x ->(scaledpT r eigenS x,mat))) $ gamma numCat alpha
                                (eigenS,mat) = quickEigen' pi s

thmmModel numCat s pi [priorZero,alpha,sigma] = ([thmm numCat pi s priors alpha sigma],[fullPi]) where
                        priors =  (replicate (numCat -1) ((1.0-priorZero)/(fromIntegral (numCat-1))) ) ++ [priorZero]
                        fullPi = makeFullPi priors pi

{-- TODO refactor this out --}
thmmModelQ numCat s pi [priorZero,alpha,sigma] = (thmmQ numCat pi s priors alpha sigma) where
                        priors = (replicate (numCat -1) ((1.0-priorZero)/(fromIntegral (numCat-1))) ) ++ [priorZero]
                        fullPi = makeFullPi priors pi


thmmPerBranchModel numCat s pi [priorZero,alpha] = ([thmmPerBranch numCat pi s priors alpha],[fullPi]) where 
                        priors = (replicate (numCat -1) ((1.0-priorZero)/(fromIntegral (numCat-1))) ) ++ [priorZero]
                        fullPi = makeFullPi priors pi
thmmPerBranchModel numCat s pi list = error $ "Fail " ++ (show list)

quickThmm numCat aln tree pi s [priorZero,alpha,sigma] = qdLkl numCat [1.0] AminoAcid aln tree (thmmModel numCat s pi) [priorZero,alpha,sigma] 
                                                        
thmm :: Int -> Vector Double -> Matrix Double -> [Double] -> Double -> Double -> BranchModel
thmm numCat pi s priors alpha sigma = (\x->(standardpT (eigQ (thmmQ numCat pi s priors alpha sigma) fullPi) x,mat)) where
                                               mat = thmmQ numCat pi s priors alpha sigma
                                               fullPi = makeFullPi priors pi

thmmQ :: Int -> Vector Double -> Matrix Double -> [Double] -> Double -> Double -> Matrix Double
thmmQ numCat pi s priors alpha sigma = fixDiag $ fromBlocks subMats where 
                                        subMats = map (\(i,mat) -> getRow i mat) $ zip [0..] qMats 
                                        qMats = (map (\(i,mat) -> setRate i mat pi) $ zip (map (*factor) $ gamma (numCat -1) alpha) $ replicate numCat $ makeQ s pi) ++ [(zeroQMat (rows s))]
                                        factor = 1.0 / (1.0 - (last priors))
                                        getRow i mat = map (kk' mat i) [0..(numCat-1)] 
                                        size = cols s
                                        kk' mat i j | i==j = mat
                                                    | otherwise = diagRect 0.0 (mapVector ((priors !! j) * sigma * ) pi) size size



thmmPerBranch :: Int -> Vector Double -> Matrix Double -> [Double] -> Double -> [Double] -> (Matrix Double,Matrix Double)
thmmPerBranch numCat pi s priors alpha [branchLength,sigma] = thmm numCat pi s priors alpha sigma [branchLength]
thmmPerBranch numCat pi s priors alpha paramList | trace ("thmmpb " ++ (show paramList)) False = undefined




loggedFunc :: (Show a) => (a -> Double) -> (a -> Double)
loggedFunc f | trace "LOGGING ON" True = f2 where
               f2 x = ans3 where
                      ans3 | trace (show x) True = ans2
                      ans2 | trace ((show x) ++ " -> " ++ (show ans)) True = ans
                      ans = f x
                    
maximize method pre maxiter size f = minimize method pre maxiter size (invert f) 

maximizeDG lower upper method pre maxiter size tol f params = minimizeD method pre maxiter size tol f2 (mygrad lower upper f2) params where
        f2 = invert f

stepSize = 1E-4

mygrad lower upper f params = mygrad' [] lower upper f params
mygrad' donep (l:lower) (u:upper) f (p:ps) = ans:(mygrad' (p:donep) lower upper f ps) where 
        f2 x = f $ (reverse (x:donep)) ++ ps
        ans = case (l,u) of 
                (Nothing,Nothing) -> fst $ derivCentral stepSize f2 p 
                (Just x,Nothing) -> fst $ derivForward stepSize f2 p
                (Nothing,Just x) -> fst $ derivBackward stepSize f2 p
                (Just li,Just ui) | (p + stepSize > ui) -> fst $ derivBackward stepSize f2 p
                (Just li,Just ui) -> fst $ derivForward stepSize f2 p
mygrad' _ _ _ _ [] = []





invert f = (*(-1)). f


optBL :: (AddTree t) => ModelF -> t -> [Double] -> [Double] -> Double -> DNode
optBL model t' params priors cutoff = traceShow optTree optTree where
        tree  = addModelFx t' (model params) priors
        optTree = optBLD $ optBLD tree

optBLx :: (AddTree t) => Int ->  ModelF -> t -> [Double] -> [Double] -> Double -> DNode
optBLx val model t' params priors cutoff = traceShow optTree optTree where
                tree  = addModelFx t' (model params) priors
                optTree = optBLD $ optBLD tree
 
genericOpt opt maxiter size precision lower upper func initialParams = fst $ maximize opt precision maxiter size (boundedFunction (-1E20) lower upper func) initialParams
simplexOpt :: [Maybe Double] -> [Maybe Double] -> ([Double] -> Double) -> [Double] -> [Double]
simplexOpt l u f p = genericOpt NMSimplex2 100000 (map (\x->0.1) p) 1e-4 l u f p
bfgsOpt l u f p = fst $ maximizeDG l u VectorBFGS2 1e-4 100000 0.1 1e-4 f p 

bobyqaOpt :: [Double] -> [Maybe Double] -> [Maybe Double] -> ([Double] -> Double) -> [Double] -> IO ([Double],Int)
bobyqaOpt ss l u f p = bobyqa ss 1E-4 p (invert f) l u
cobylaOpt ss l u f p = cobyla ss 1E-4 p (invert f) l u


getFunc1 priors model tree params = logLikelihood $ getFuncT1 priors model tree params

-- only parameters
getFuncT1 :: [Double] -> ModelF -> DNode -> [Double] -> DNode
getFuncT1 priors model tree = func where
        func params = addModelFx tree (model params) priors

-- bsParams then parameters 

getFuncT1A  :: [Double] -> (DNode -> Int) -> (Int,Int) -> ModelF -> DNode -> [Double] -> DNode
getFuncT1A priors mapping (0,0) model tree = getFuncT1 priors model tree
getFuncT1A priors mapping numBSParam model tree = func where
        func params = addModelFx t2 (model params') priors where
                (paramsPerBranch,paramCats) = numBSParam
                t2 = setBLMapped 1 tree getParamF 
                perBranchParams = take paramCats $ splitLists paramsPerBranch params
                params' = drop (paramsPerBranch * paramCats) params
                getParamF node | trace ("params' " ++ (show params')) True = (perBranchParams !! (mapping node))

getParamsT1 params tree = params

-- branchLengths, then other parameters
getFuncT2 :: [Double] -> ModelF -> DNode -> ([Double] -> DNode)
getFuncT2 priors model tree | trace "OK1" True = func where
        func params | trace "OK2" True  = addModelFx t2 (model params2) priors where
            (t2,params2') = setBLX' 0 (map (\x->[x]) params) tree
            params2 = map head params2'

getParamsT2 params tree = (map head $ getBL' tree []) ++  params

-- branchLengths, then bsParams, then other params
getFuncT3 :: [Double] -> (DNode -> Int) -> (Int,Int) -> ModelF -> DNode -> ([Double] -> DNode)
getFuncT3 priors mapping (0,0) model tree = getFuncT2 priors model tree
getFuncT3 priors mapping (paramsPerBranch,paramCats) model tree = func where
        func params = addModelFx t2 (model params2) priors where
          indicies  = getLinearMap mapping tree []
          bsParam = splitLists paramsPerBranch $ drop (length indicies) params
          branchLengths = map (\(x,y) -> y:(bsParam !! x)) $ zip indicies params 
          (t2,_)  = setBL' branchLengths tree
          params2 = drop ((length indicies) + paramsPerBranch * paramCats) params

getParamsT3 = getParamsT2 

splitLists i [] = []
splitLists i x@(z:zs) = y:(splitLists i ys) where (y,ys) = splitAt i x
                 
data IterType = Full | BL | Params deriving (Eq, Show)

optBSParamsBL numBSParam mapping = optWithBS' [] 1E-4 numBSParam (Just mapping)
optBSParamsBLIO numBSParam mapping initialStepSize = optWithBSIO' [] 1E-4 numBSParam (Just mapping) initialStepSize (map (/100.0) initialStepSize)
optParamsBL = optWithBS' [] 1E-4 (0,0) Nothing


optWithBSIO' :: [(Double,IterType)] -> Double -> (Int,Int) -> (Maybe (DNode -> Int)) -> [Double] -> [Double] -> [Maybe Double] -> [Maybe Double] -> [Double] -> ModelF -> DNode -> [Double] -> IO (DNode,[Double])
optWithBSIO' iterations cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree startParams = do 
     ans <- case iterations of
                ((x,t):(x',t'):(x'',t''):(x''',t'''):xs) | (x-x''' < cutoff) && (t'''==Full || t==Full) && stepSizeMet -> return (tree,startParams) --stop
                ((x,t):(x',t'):(x'',t''):(x''',t'''):xs) | (x-x''' < cutoff) && (t'''==Full || t==Full) -> optWithBSIO' iterations cutoff numBSParam mapping (incrementStepSize stepSize) limitStepSize lower upper priors model tree startParams
                list@((x,BL):xs) -> do 
                           (bestParams',err) <- bobyqaOpt stepSize lower upper (loggedFunc $ logLikelihood . getFuncT1A priors (fromJust mapping) numBSParam model tree) (enforceBounds lower upper startParams)
                           let tree' = getFuncT1A priors (fromJust mapping) numBSParam model tree bestParams'
                           let lkl = logLikelihood tree'
                           optWithBSIO' ((lkl,Params):list) cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree' bestParams' where -- Opt Params
                list@((x,t):(x',t'):_:_:xs) | (x-x' < (cutoff*1000)) && (t /= Full) -> do 
                           let myBL = map head $ getBL' tree [] 
                           let numBL = length myBL
                           let addBL p = myBL ++ startParams
                           let dropBL p = drop numBL p
                           let lower' = (replicate numBL $ Just 1E-6) ++ lower
                           let upper' = (replicate numBL $ Just 10.0) ++ upper
                           let stepSize' = (replicate numBL $ 0.01) ++ (map (/2) stepSize)
                           {-let startParams' = enforceBounds lower' upper' (addBL startParams)-}
                           let startParams' = addBL startParams
                           (bestParams',err) <- bobyqaOpt stepSize' lower' upper' (loggedFunc $ logLikelihood . getFuncT3 priors (fromJust mapping) numBSParam model tree) startParams'
                           --let (bestParams',err) = output
                           let tree' = getFuncT3 priors (fromJust mapping) numBSParam model tree bestParams'
                           let lkl = logLikelihood tree'
                           optWithBSIO' ((lkl,Full):list) cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree' (dropBL bestParams') where --full Opt
                list -> trace "BL OPT" $ optWithBSIO' ((lkl,BL):list) cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree' startParams where -- Opt Params
                                tree' = optBLDFull0 $ getFuncT1A priors (fromJust mapping) numBSParam model tree startParams
                                params' = drop ((fst numBSParam) * (snd numBSParam)) startParams
                                lkl = logLikelihood $ getFuncT1A priors (fromJust mapping) numBSParam model tree' startParams
     return ans where
     incrementStepSize ss = map (\(x,y) -> max x y) $ zip limitStepSize $ map (/2) ss
     stepSizeMet = stepSizeMet' limitStepSize stepSize
     stepSizeMet' a b = (Nothing /=) $ find (\(x,y)-> (<) x y) $ zip a b
     

 
optWithBS' :: [(Double,IterType)] -> Double -> (Int,Int) -> (Maybe (DNode -> Int)) -> [Maybe Double] -> [Maybe Double] -> [Double] -> ModelF -> DNode -> [Double] -> (DNode,[Double])
optWithBS' iterations cutoff numBSParam mapping lower upper priors model tree startParams = trace ((show (take 10 iterations)) ++ "\n" ++ (show startParams) ++ "\n" ++ (show tree)) (bestTree,bestParams) where
        (bestTree,bestParams) = case iterations of
                ((x,t):(x',t'):(x'',t''):xs) | (x-x'' < cutoff) && (t''==Full || t==Full) -> (tree,startParams) --stop
                list@((x,t):(x',t'):_:_:xs) | (x-x' < 0.5) && (t /= Full) -> trace "FULLOPT" $  optWithBS' ((lkl,Full):list) cutoff numBSParam mapping lower upper priors model tree (dropBL bestParams') where --full Opt
                            bestParams' = simplexOpt lower' upper' (loggedFunc $ logLikelihood . getFuncT3 priors (fromJust mapping) numBSParam model tree) startParams'
                            startParams' = enforceBounds lower' upper' (addBL startParams)
                            tree' = getFuncT3 priors (fromJust mapping) numBSParam model tree bestParams'
                            lkl = logLikelihood tree'
                            myBL = map head $ getBL' tree [] 
                            numBL = length myBL
                            addBL p = myBL ++ startParams
                            dropBL p = drop numBL p
                            lower' = (replicate numBL $ Just 0.001) ++ lower
                            upper' = (replicate numBL $ Just 10.0) ++ upper
                list@((x,BL):xs) -> trace "PARAM OPT" $ optWithBS' ((lkl,Params):list) cutoff numBSParam mapping lower upper priors model tree' bestParams' where -- Opt Params
                                bestParams' | trace "OK" True = simplexOpt lower upper (loggedFunc $ logLikelihood . getFuncT1A priors (fromJust mapping) numBSParam model tree) (enforceBounds lower upper startParams)
                                lkl = logLikelihood tree'
                                tree' = getFuncT1A priors (fromJust mapping) numBSParam model tree bestParams'
                list -> trace "BL OPT" $ optWithBS' ((lkl,BL):list) cutoff numBSParam mapping lower upper priors model tree' startParams where -- Opt Params
                                tree' = optBLDFull0 $ getFuncT1A priors (fromJust mapping) numBSParam model tree startParams
                                params' = drop ((fst numBSParam) * (snd numBSParam)) startParams
                                lkl = logLikelihood $ getFuncT1A priors (fromJust mapping) numBSParam model tree' startParams
                

enforceBounds l u p = enforce (>) u $ enforce (<) l p where
        enforce f (x:xs) (y:ys) = (enforce' f x y):(enforce f xs ys)
        enforce f [] [] = []
        enforce' f x y = case x of
                        Nothing -> y
                        Just x' -> if (f y x') then x' else y


--TODO eradicate this function!
dummyTree :: (AddTree t) => t -> DNode 
dummyTree t = addModelFx t (basicModel wagS wagPi []) [1.0]

whichModels t mapping = map mapping (nonRootNodes $ dummyTree t)

makeMapping :: (AddTree t) => (([String],[String]) -> a) -> t -> (DNode -> a) 
makeMapping splitmap t = lookupF where
        splits = getSplits $ dummyTree t
        splitList = zip splits $ map splitmap splits 
        eitherSideList = splitList >>= (\((a,b),c)->[(a,c),(b,c)]) 
        lookupF node = case find (\(x,a) -> (sort x) == (sort $ descendents node)) eitherSideList of
                            Just (_,y) -> y
                            Nothing -> error $ "Can't find in " ++ (show (sort $ descendents node)) ++ "\nin " ++ (show $ map fst eitherSideList) 


traceX x = traceShow x x
traceXP p x = trace (p ++ " " ++ (show x)) x


genList :: StdGen -> [StdGen]
genList stdGen = first:remainder where
        (first,next) = split stdGen
        remainder = genList next

makeSimulatedTree seqDataType hiddenClasses stdGen length (DTree l m r _ _ priors pis) = DTree left middle right pLs (replicate length 1) priors pis where
        left = makeSimulatedTree' seqDataType hiddenClasses leftR topSeq models l
        middle = makeSimulatedTree' seqDataType hiddenClasses middleR topSeq models m
        right = makeSimulatedTree' seqDataType hiddenClasses rightR topSeq models r
        pLs = calcRootPL l m r
        (r0:r1:leftR:middleR:rightR:remainder) = genList stdGen
        models = take length $ map (drawFromDist priors) $ randomRs (0.0,1.0) r0
        drawLetters :: [Int]
        drawLetters = map (\(p,pri) -> drawFromDist pri p) $ zip (randomRs (0.0,1.0) r1) $ map toList $ map (pis!!) models
        topSeq :: [Int]
        topSeq = take length drawLetters 

makeSimulatedTree' :: SeqDataType -> Int -> StdGen -> [Int] -> [Int] -> DNode -> DNode
makeSimulatedTree' Nucleotide _ _ _ _ _ = error "Nucleotide simulation unimplemented"

makeSimulatedTree' AminoAcid hiddenClasses stdGen topSeq models (DLeaf name dist _ _ modelList _) = DLeaf name dist bottomSeq tipLkl modelList pLs where
        myMats :: [Matrix Double]
        myMats = map (\x -> fst $x dist) modelList
        myVectors :: [[[Double]]]
        myVectors = map toLists $ map (myMats!!) models
        myVectors' = (map (\(vec,base) -> vec!!base) $ zip myVectors topSeq)
        bottomSeq = map (aaOrder!!) $ map (`mod` 20 ) $ map (\(p,pris) -> drawFromDist pris p) $ zip (randomRs (0.0,1.0) stdGen) myVectors'
        pLs = calcLeafPL tipLkl dist modelList
        tipLkl = getPartial hiddenClasses AminoAcid bottomSeq

makeSimulatedTree' seqDataType hiddenClasses stdGen topSeq models (DINode l r dist modelList _)  = DINode left right dist modelList pLs where
        myMats :: [Matrix Double]
        myMats = map (\x -> fst$ x dist) modelList
        myVectors :: [[[Double]]]
        myVectors = map toLists $ map (myMats!!) models
        (r0:leftR:rightR:_) = genList stdGen 
        bottomSeqPartial :: [(Double,[Double])]
        bottomSeqPartial = zip (randomRs (0.0,1.0) r0) (map (\(vec,base) -> vec!!base) $ zip myVectors topSeq) 
        bottomSeq = map (\(p,pris) -> drawFromDist pris p) bottomSeqPartial
        left = makeSimulatedTree' seqDataType hiddenClasses leftR bottomSeq models l
        right = makeSimulatedTree' seqDataType hiddenClasses rightR bottomSeq models r
        pLs = calcPL left right dist modelList

drawFromDist pris = drawFromDist' $ map nonNeg pris where
        nonNeg x | x < 0.0 = 0.0
        nonNeg x = x

drawFromDist' pris | (sum pris) < 1.0 = drawFromDist' normPris where
        normPris = map (/total) pris
        total = sum pris

drawFromDist' pris = drawFromDist'' 0 pris

drawFromDist'' :: Int -> [Double] -> Double -> Int
--drawFromDist' c pris p | trace ((show c) ++ " " ++ (show pris) ++ " " ++ (show p)) False = undefined
drawFromDist'' c (pri:pris) p = if (p <= pri) then c else (drawFromDist'' (c+1) pris (p-pri))

simulateSequences seqDataType hiddenClasses stdGen length root = quickListAlignment names seqs where
        simTree  = makeSimulatedTree seqDataType hiddenClasses stdGen length root
        leaves = getLeaves simTree
        names = map dName leaves
        seqs = map Phylo.Likelihood.sequence leaves


annotateTreeWithNumberSwitchesSigma sigma (DTree l m r models patcounts priors pis) = DTree (an l) (an m) (an r) models patcounts priors pis where
                                                                            an = annotateTreeWithNumberSwitchesSigma' sigma priors pis

annotateTreeWithNumberSwitchesSigma' sigma priors pis (DINode l r (bl:[]) models mats) = DINode (an l) (an r) [bl,sigma,switchbl] models mats where
                                                                                        switchbl = calcSwitchBL priors pis models [bl,sigma]
                                                                                        an = annotateTreeWithNumberSwitchesSigma' sigma priors pis

annotateTreeWithNumberSwitchesSigma' sigma priors pis (DLeaf name (bl:[]) seq tip models partial) = DLeaf name [bl,sigma,switchbl] seq tip models partial where
                                                                                        switchbl = calcSwitchBL priors pis models [bl,sigma]

annotateTreeWithNumberSwitches (DTree l m r models patcounts priors pis) = DTree (an l) (an m) (an r) models patcounts priors pis where
                                                                            an = annotateTreeWithNumberSwitches' priors pis
annotateTreeWithNumberSwitches' priors pis (DINode l r (bl:sigma:[]) models mats) = DINode (an l) (an r) [bl,sigma,switchbl] models mats where
                                                                                        switchbl = calcSwitchBL priors pis models [bl,sigma]
                                                                                        an = annotateTreeWithNumberSwitches' priors pis

annotateTreeWithNumberSwitches' priors pis (DLeaf name (bl:sigma:[]) seq tip models partial) = DLeaf name [bl,sigma,switchbl] seq tip models partial where
                                                                                        switchbl = calcSwitchBL priors pis models [bl,sigma]


calcSwitchBL priors pis models (bl:sigma:[]) = (*) bl $ sum $ zipWith (*) priors switchingrates where
                                                 mats = map (\x-> snd $ x [bl,sigma]) models
                                                 switchingrates = map (switchingSum 20) $ zip pis mats --fix this magic number

switchingSum nc (pi,mat) = getSwitchingRate mat pi nc

                                                                                

