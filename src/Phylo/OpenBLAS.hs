{-# LANGUAGE ForeignFunctionInterface #-}
module Phylo.OpenBLAS where
import Foreign
import Foreign.C.Types

foreign import ccall safe "goto_set_num_threads" 
        c_set_num_threads :: CInt -> IO ()

set_num_threads :: Int -> IO ()
set_num_threads i = c_set_num_threads (fromIntegral i)
