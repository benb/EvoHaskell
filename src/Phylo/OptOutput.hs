module Phylo.OptOutput where
import Phylo.Likelihood
import Control.DeepSeq
import Text.Printf
import System.IO (hFlush,stdout)


printNiceOpt len (x,x') = do
   let string = getNiceOpt2 len x x'
   let cr = string `deepseq` '\r'
   putChar cr
   putStr (replicate len ' ')
   putChar cr
   putStr string
   hFlush stdout

getNiceOpt len (a,b,c) = outL ++ outR where
                                outL = (show c) ++ " : " ++ (format $ logLikelihood a)
                                outR = replicate (len - (length outL)) ' '
                                format = printf "%.5f"

getNiceOpt2 len (a,b,c) (a',b',c') = outL ++ outR where
                                        outL = (show c') ++ " : " ++ (format $ logLikelihood a) ++ " -> " ++ (format $ logLikelihood a')
                                        format = printf "%.5f"
                                        outR = replicate (len - (length outL)) ' '


optWithOutput ans = do putStrLn "optimising model"
                       let start = head ans
                       putStr $ getNiceOpt 100 $ head ans
                       hFlush stdout
                       mapM (\x -> do printNiceOpt 100 x) (zip ans (tail ans))
                       let (a,b,_) = last ans
                       putStrLn ""
                       return (cachedBranchModelTree a,b)

