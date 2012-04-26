module Phylo.Graphics.Plotting where
import Graphics.Rendering.Chart
import Data.Accessor                                                                                                                                                                                                                          
import Data.Accessor.Tuple                                                                                                                                                                                                                    
import Data.Colour                                                                                                                                                                                                                            
import Data.Colour.Names                                                                                                                                                                                                                      
import Data.Colour.SRGB 
import Debug.Trace

data OutputType = PNG | PS | PDF | SVG                                                                                                                                                                                                        

chooseLineWidth PNG = 1.0                                                                                                                                                                                                                     
chooseLineWidth PDF = 0.25                                                                                                                                                                                                                    
chooseLineWidth PS = 0.25                                                                                                                                                                                                                     
chooseLineWidth SVG = 0.25             

makePlot :: [Double] -> [(Double,Double)] -> OutputType -> Renderable ()
makePlot qD nf otype = toRenderable layout where                                                                                                                                                                                                 
    layout = layout1_title ^=""                                                                                                                                                                                                         
           $ layout1_background ^= solidFillStyle (opaque white)                                                                                                                                                                              
           $ layout1_left_axis ^: laxis_generate ^= const baxis
           $ layout1_bottom_axis ^: laxis_title ^= "theoretical"
           $ layout1_left_axis ^: laxis_title ^= "empirical"
           $ layout1_bottom_axis ^: laxis_generate ^= const baxis
           $ layout1_plots ^= [ Left (toPlot xEqualYLine)
                              , Left (toPlot qDPoints)                                                                                                                                                                                          
                              , Left (toPlot qDArea)                                                                                                                                                                                       
                              ]                                                                                                                                                                                                               
           $ setLayout1Foreground fg                                                                                                                                                                                                          
           $ defaultLayout1                                                                                                                                                                                                                   

    fg = opaque black 
    qLen = length qD
    xvals =  [(fromIntegral x)/(fromIntegral qLen+1)|x<-[1..(qLen+1)]]
    xyvals :: [(Double,Double)]
    xyvals = zip xvals qD
    xyArea = zip xvals nf
    xEqualY = zip xvals xvals

    xEqualYLine = plot_lines_style ^= dashedLine 1.0 [1.0,1.0] (opaque black)
                $ plot_lines_values ^= [xEqualY]
                $ defaultPlotLines

    qDLine = plot_lines_style ^= lineStyle (opaque blue)                                                                                                                                                                                      
           $ plot_lines_values ^= [xyvals]                                                                                                                                                                                            
           $ plot_lines_title ^= ""                                                                                                                                                                                                  
           $ defaultPlotLines            

    qDPoints = plot_points_style ^= filledCircles 2 (opaque blue)                                                                                                                                                                                
           $ plot_points_values ^= xyvals
           $ plot_points_title ^= ""                                                                                                                                                                                                 
           $ defaultPlotPoints   

    qDArea = plot_fillbetween_values ^= xyArea                                                                                                                                                 
                $ plot_fillbetween_style  ^= solidFillStyle (withOpacity blue 0.2)                                                                                                                                                            
                $ defaultPlotFillBetween    

    lineStyle c = line_width ^= 3 * chooseLineWidth otype                                                                                                                                                                                     
                $ line_color ^= c                                                                                                                                                                                                             
                $ defaultPlotLines ^. plot_lines_style 
    baxis :: AxisData Double

    baxis = AxisData {                                                                                                                                                                                                                        
        axis_viewport_ = vmap (0.0,1.0),
        axis_tropweiv_ = invmap (0.0,1.0),
        axis_ticks_    = [(v/10,3) | v <- [0,2..10]],
        axis_grid_     = [0.0,0.1..1.0],
        axis_labels_   = [[(v/10,show (v/10)) | v <- [0,2..10]]]
    }
