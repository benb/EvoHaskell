{-# LANGUAGE ForeignFunctionInterface #-}
module Phylo.NLOpt (bobyqa,cobyla) where

import Data.Maybe (fromMaybe)
import Foreign
import Foreign.C
import Foreign.C.Types
import Foreign.Ptr (Ptr, FunPtr, freeHaskellFunPtr)
import Control.Monad
import Debug.Trace
import Control.Applicative
import Data.Maybe

foreign import ccall safe "nlopt_c opt_bobyqa" 
        bobyqa_ :: CDouble -> Ptr CDouble -> Ptr CDouble -> CUInt -> FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)) -> Ptr CDouble -> Ptr CDouble -> CInt
foreign import ccall safe "nlopt_c opt_cobyla" 
        cobyla_ :: CDouble -> Ptr CDouble -> Ptr CDouble -> CUInt -> FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)) -> Ptr CDouble -> Ptr CDouble -> CInt
foreign import ccall safe "nlopt_c opt_mma"
        mma_ :: CDouble -> Ptr CDouble -> Ptr CDouble -> CUInt -> FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)) -> Ptr CDouble -> Ptr CDouble -> CInt

foreign import ccall "wrapper"
        wrap :: (CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble) -> IO (FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)))

nlopt :: (CDouble -> Ptr CDouble -> Ptr CDouble -> CUInt -> FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)) -> Ptr CDouble -> Ptr CDouble -> CInt) -> [Double] -> Double ->  [Double] -> ([Double] -> (Double,Maybe [Double])) -> [Maybe Double] -> [Maybe Double] -> IO ([Double],Int)

bobyqa = nlopt bobyqa_
cobyla = nlopt cobyla_
traceX a x = trace (show a ++ (show x)) x
traceXm a = liftM (traceX a)
nlopt met stepSize xtol params f lower upper = do lower' <- newArray $ map (realToFrac . fromMaybe (-1E100)) lower
                                                  upper' <- newArray $ map (realToFrac . fromMaybe 1E100) upper 
                                                  stepSize' <- newArray $ map realToFrac stepSize
                                                  let np = trace ("Start params " ++ (show params) ++ " np " ++ (show (length params))) $ length params
                                                  let f' a b c d = do (ans,deriv)<-(liftM f) (fmap (map realToFrac) $ peekArray np b)
                                                                      case nullPtr of 
                                                                            x | x==nullPtr -> return ()
                                                                            ptr -> pokeArray ptr $ fromJust deriv
                                                                      return $ realToFrac ans
                                                  f'' <- wrap f'
                                                  startP <- newArray $ traceX ("realToFrac") $ map (realToFrac) params
                                                  let retCode = fromIntegral $ met (realToFrac xtol) stepSize' startP (fromIntegral np) f'' lower' upper' 
                                                  ans <- seq retCode $ fmap (map realToFrac) $ peekArray np startP
                                                  freeHaskellFunPtr f''
                                                  return (ans,retCode)

