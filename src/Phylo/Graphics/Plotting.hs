module Phylo.Graphics.Plotting where
import Graphics.Rendering.Chart
import Data.Accessor                                                                                                                                                                                                                          
import Data.Accessor.Tuple                                                                                                                                                                                                                    
import Data.Colour                                                                                                                                                                                                                            
import Data.Colour.Names                                                                                                                                                                                                                      
import Data.Colour.SRGB 

data OutputType = PNG | PS | PDF | SVG                                                                                                                                                                                                        

chooseLineWidth PNG = 1.0                                                                                                                                                                                                                     
chooseLineWidth PDF = 0.25                                                                                                                                                                                                                    
chooseLineWidth PS = 0.25                                                                                                                                                                                                                     
chooseLineWidth SVG = 0.25             

makePlot :: [Double] -> OutputType -> Renderable ()
makePlot qD otype = toRenderable layout where                                                                                                                                                                                                 
    layout = layout1_title ^="QQPlot"                                                                                                                                                                                                         
           $ layout1_background ^= solidFillStyle (opaque white)                                                                                                                                                                              
           $ layout1_left_axis ^: laxis_generate ^= const baxis
           $ layout1_bottom_axis ^: laxis_generate ^= const baxis
           $ layout1_plots ^= [ Left (toPlot qDPoints)                                                                                                                                                                                          
                              , Left (toPlot qDArea)                                                                                                                                                                                       
                              ]                                                                                                                                                                                                               
           $ setLayout1Foreground fg                                                                                                                                                                                                          
           $ defaultLayout1                                                                                                                                                                                                                   

    fg = opaque black 
    xvals = [0.0,0.01..]
    xyvals :: [(Double,Double)]
    xyvals = zip xvals qD

    qDLine = plot_lines_style ^= lineStyle (opaque blue)                                                                                                                                                                                      
           $ plot_lines_values ^= [xyvals]                                                                                                                                                                                            
           $ plot_lines_title ^= "empirical"                                                                                                                                                                                                  
           $ defaultPlotLines            

    qDPoints = plot_points_style ^= filledCircles 2 (opaque blue)                                                                                                                                                                                
           $ plot_points_values ^= xyvals
           $ plot_points_title ^= "empirical"                                                                                                                                                                                                 
           $ defaultPlotPoints   

    qDArea = plot_fillbetween_values ^= [(d, (v * 0.95, v * 1.05)) | (d,v) <-xyvals]                                                                                                                                                  
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
