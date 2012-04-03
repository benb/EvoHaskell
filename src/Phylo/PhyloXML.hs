module Phylo.PhyloXML where
import Phylo.Likelihood 
import Debug.Trace
import Data.List
import Data.Colour.SRGB
import Data.Colour.Names
import Data.Colour

data BranchLengthType = BranchLength | Custom String

branchLengthSigma = Custom "Sigma"
branchLengthSwitching = Custom "Switching"


xmlHeader = ["<?xml version=\"1.0\" encoding=\"UTF-8\"?>","<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\">"]
phyloXML taginfo tree = xmlHeader ++ (xmlTree taginfo tree) ++ ["</phyloxml>"]
simplePhyloXML = phyloXML [Just BranchLength, Nothing, Just branchLengthSwitching] 
quantilePhyloXML = phyloXML [Just BranchLength, Just $ Custom "empirical_quantile"]
xmlTree ti (DTree l m r _ _ _ _ ) = ["<phylogeny rooted=\"false\">"] ++ ["<clade>"] ++ (xmlTree ti l) ++ (xmlTree ti m) ++ (xmlTree ti r) ++ ["</clade>"] ++ ["</phylogeny>"]
xmlTree ti (DINode l r bl _ _ ) = ("<clade>":(calcBL ti bl)) ++ (xmlTree ti l) ++ (xmlTree ti r) ++ ["</clade>"]
xmlTree ti (DLeaf name bl _ _ _ _ ) = ["<clade>","<name>"++name++"</name>"] ++ (calcBL ti bl) ++ ["</clade>"]

calcBL taginfo bl = map snd $ sortBy order $ concatMap xmlise $ zip taginfo bl where
                        xmlise a = case a of
                                 (Nothing,x) -> []
                                 (Just BranchLength,x) -> [(1,"<branch_length>" ++ (show x) ++ "</branch_length>")]
                                 (Just (Custom name),x) -> [(13,"<property applies_to=\"parent_branch\" ref=\"man:" ++ name ++ "\" datatype=\"xsd:double\">" ++ (show x) ++ "</property>")
                                                           --,(4,"<width>" ++ (show x) ++ "</width>")
                                                           ,(5,colorise x)
                                                           ]
                        order a b = compare (fst a) (fst b)

colorise x = let blended = toSRGB24 $ blend x red blue in 
                 intercalate "\n" ["<color>","<red>" ++ (show $ channelRed blended) ++ "</red>"
                                  ,"<green>" ++ (show $channelGreen blended) ++ "</green>"
                                  ,"<blue>" ++ (show $ channelBlue blended) ++ "</blue>"
                                  ,"</color>"]
