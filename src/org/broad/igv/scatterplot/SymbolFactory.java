package org.broad.igv.scatterplot;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;


/*
*   Factory class for constructing plot symbol settings for series plots.
*
*   Note: scatter Plot data container class is necessary to determine the number
*   of series to set up for, when al are set up at once.
*
* @author martind
*
* */

 public class SymbolFactory {

    // plot shape types - Ring and Oval always have outlines
    // SeriesSpecific is for a different shape per series, in the order defined here
    public enum IGVSymbolShape { SeriesSpecific, Circle, Triangle, Rectangle, Ellipse,
          Diamond, Gradient, Pentagon, Hexagon, Octagon, Ring, Oval,  Square }

    // default colors assigned per series, starting with systemColor[1] (red);
    // systemColors[0] (black) is reserved for outlining shapes
     private Color [] defaultColor = {Color.black, Color.red,  Color.green, Color.blue, Color.yellow,
        Color.magenta, Color.cyan, Color.orange, Color.pink, Color.gray, Color.darkGray,
        Color.lightGray, Color.white};

    // default plot shape is set to a Rectangel with width = 2, height = 2
    private IGVSymbolShape igvDefaultSymbolShape;   // default symbol shape
    private int igvShapeWidth;          // default shape witdh
    private int igvShapeHeight;         // default shape height
    private boolean igvShapeOutlined;   // shape has outline; false for no outline
    private boolean igvShapeFilled;     // shape is filled; false for no fill

    // IGV data container class indicates attributes and categories per attribute
    //private ArrayList<String> igvSymbolCategories; // symbol categories for series labeling
    private ArrayList<Shape> igvSymbolShapeList;        // symbol shapes for series plots
    private ArrayList<Color> igvSymbolColorList;        // symbol colors for series plots

    // symbol settings constructed for either default or supplied shape specifications
    private ArrayList<SymbolSettings> igvSymbolSettingsList;

    /*
    *   Null constructor is provided for caller defined symbol specifications.
    *
    *   Note: Symbol specifications are provided through the addSymbol method.
    */
    public SymbolFactory() {

        // set up empty symbol settings array for user specified
        igvSymbolSettingsList = new ArrayList<SymbolSettings>();
    }

    /*
    *   Constructor sets up symbol settings with default shapes and colors.
    *
    *   parameters:
    *       categories - attribute categories to be assigned symbol settings
    *       shapeWidth -  shape bounding rectangle longer dimension
    *       shapeHeight -  shape bounding rectangle height
    *       isOutlined - boolean flag is true for outlined shape; false for none
    *       isFilled - boolean flag is true for filled shape; false for hollow
    * */
    public SymbolFactory(String[] categories, int shapeWidth,
                         int shapeHeight, boolean isOutlined, boolean isFilled){

        igvDefaultSymbolShape = IGVSymbolShape.SeriesSpecific;
        igvShapeWidth = shapeWidth;
        igvShapeHeight = shapeHeight;
        igvShapeOutlined = isOutlined;
        igvShapeFilled = isFilled;

        createDefaultSymbolSettings(categories);
    }

    /*
    *   Constructor sets up symbol factory with default symbol information.
    *
    *   parameters:
    *       categories - attribute categories to be assigned symbol settings.
    *       defaultShape - shape selected for default series shapes
    *       shapeWidth -  shape bounding rectangle longer dimension
    *       shapeHeight -  shape bounding rectangle height
    *       isOutlined - boolean flag is true for outlined shape; false for none
    *       isFilled - boolean flag is true for filled shape; false for hollow
    * */
    public SymbolFactory(String[] categories, IGVSymbolShape defaultShape, int shapeWidth,
                         int shapeHeight, boolean isOutlined, boolean isFilled){

        igvDefaultSymbolShape = defaultShape;
        igvShapeWidth = shapeWidth;
        igvShapeHeight = shapeHeight;
        igvShapeOutlined = isOutlined;
        igvShapeFilled = isFilled;

        createDefaultSymbolSettings(categories);
    }

    public IGVSymbolShape getDefaultSymbolShape() {
        return igvDefaultSymbolShape;
    }

    public int getDefaultSymbolWidth() {
        return igvShapeWidth;
    }

    public int getIgvShapeHeight() {
        return igvShapeHeight;
    }

    public boolean isDefaultShapeOutlined() {
        return igvShapeOutlined;
    }

    public boolean isDefaultShapeFilled() {
        return igvShapeFilled;
    }
    /*
    *   Method exports the series symbol settings constructed by the factory.
    *
    *   Returns:
    *       ArrayList of SymbolSettings which define the plotting symboles for
    *       the categories the factory was constructed with, or had added to it.
    *
    * */
    public SymbolSettings[] getSeriesSymbolSettings(){
        int nSeries = igvSymbolSettingsList.size();

        // test if series symbols were loaded
        if(igvSymbolSettingsList.size() == 0)
            return null;
        
        SymbolSettings[] seriesSymbolSettings = new SymbolSettings[nSeries];
        
        return igvSymbolSettingsList.toArray(seriesSymbolSettings);
    }

    /*
    *   Method adds caller defined shapes to the symbol settings ArrayList
    *   for series plot settings.
    *
    *   Note: method is static so it can be used as a simple factory tool
    *       to create symbol settings for plotting a particular shape.
    *
    *   Parameters:
    *       series - series number requiring a shape (starting index is 0)
    *       shape - shape for the series plot
    *       color - shape color for fill or outline (filled outline is black)
    *       isOutlined - boolean indicates if shape has an outline
    *       isFilled - boolean indicates shape is filled or unfilled
    *
    *   Note: Throws IndexOutOfBoundsException if series > current size of the symbol Settings;
    *       i.e. series entry is either already defined, or one greater than the last previously defined.
    *
    * */
    public void addSeriesSymbol(int series, String category, Shape shape, Color color,
                                boolean isOutlined, boolean isFilled){

         boolean isVisible = true; // assuming caller wants all series initially visible

         SymbolSettings symbolSettings = new SymbolSettings(series, category,
                 igvSymbolShapeList.get(series), igvSymbolColorList.get(series), isOutlined,
                 isFilled, isVisible);

         igvSymbolSettingsList.add(series, symbolSettings);
    }

    /*
    *   Method generates a symbol shape for series plots, and returns it.
    *
    *   Note: Method is static so it can be used alone as a tool
    *       to create plotting shapes.
    *
    *   Parameters:
    *        symbolShape - shape selection for the series plot
    *        width - width of the shape's bounding rectangle
    *        height - height of the shape's bounding rectangle
    *        outlineThickness - width for outlined shapes
    *        isFilled - boolean indicates shape is filled or unfilled
    *
    *   retuns:
    *       Shapes contains all shape plot information needed for rendering
    *
    * */
    static Shape createSymbolShape(IGVSymbolShape symbolShape, int width, int height){
        Shape shape;

        // Note: xPosition, yPosition is top left of shape bounding rectangle
        // and plot location for scatterplot is mapped to the center (0,0)
        float xPosition;
        float yPosition;

        if(symbolShape == IGVSymbolShape.Circle) {
            float radius = Math.min(width, height);
            xPosition = -radius/2.0f;
            yPosition = -radius/2.0f;
            shape = new Arc2D.Float(xPosition, yPosition, (float) width, (float) height,
                    0.0f, 360.0f, Arc2D.OPEN);
        }

        else if(symbolShape == IGVSymbolShape.Ellipse) {
            xPosition = -width /2.0f;
            yPosition = -height /2.0f;
            shape = new Ellipse2D.Float(xPosition, yPosition, (float) width, (float) height);
        }

        else if(symbolShape == IGVSymbolShape.Triangle) {
            int baseX = (int)(Math.min(width, height)/2.0f);
            int baseY = (int)(Math.max(width, height)/2.0f);
            int[] xPoints = {-baseX, 0, baseX};
            int[] yPoints = {baseY, -baseY, baseY};
            shape = new Polygon(xPoints, yPoints, 3);
        }

        else if(symbolShape == IGVSymbolShape.Square) {
            width = Math.min(width, height);
            xPosition = -width /2.0f;
            yPosition = -height /2.0f;
            shape = new Rectangle2D.Float(xPosition, yPosition, width, height);
        }

        else if(symbolShape == IGVSymbolShape.Rectangle) {
            xPosition = -width /2.0f;
            yPosition = -height /2.0f;
            shape = new Rectangle2D.Float(xPosition, yPosition, width, height);
        }

        else if(symbolShape == IGVSymbolShape.Diamond){
            int baseX = (int)(Math.min(width, height)/2.0f);
            int baseY = (int)(Math.max(width, height)/2.0f);
            int[] xPoints = {-baseX, 0, baseX, 0};
            int[] yPoints = {0, -baseY, 0, baseY};
            shape = new Polygon(xPoints, yPoints, 4);
        }

        else if(symbolShape == IGVSymbolShape.Gradient){
             int baseX = (int)(Math.min(width, height)/2.0f);
            int baseY = (int)(Math.max(width, height)/2.0f);
            int[] xPoints = {-baseX, 0, baseX};
            int[] yPoints = {-baseY, baseY, -baseY};
            shape = new Polygon(xPoints, yPoints, 3);
        }

        // symbol settings for squares and rectangle
        // Note: Need a shape at least 4 x 4
        else if(symbolShape == IGVSymbolShape.Octagon) {
            int baseX = (int)(Math.max(width, height) /4.0f);
            int baseY = (int)(Math.max(width, height)/4.0f);
            int[] xPoints = {-baseX * 2, -baseX, baseX,  baseX * 2,
                baseX * 2, baseX, -baseX, -baseX * 2};
            int[] yPoints = {-baseY, -baseY * 2, -baseY * 2, -baseY,
                baseY, baseY * 2, baseY * 2, baseY};
            shape = new Polygon(xPoints, yPoints, 8);
        }

        // symbol not defined
        else
            shape = null;

        return shape;
    }

    /*
    *   Factory method for creating attribute symbol settings with the
    *   default shapes, colors, and display attributes.
    *
    *   Parameters:
    *       categories - category names for plot series.
    *
    *   Sets up igvSymbolSettingsList ArrayList
    *
    * */
    private void createDefaultSymbolSettings(String[] categories){

        SymbolSettings symbolSettings;

       // Number of series will be number of categories for an attribute.
       int nSeries = categories.length;
       setupSymbolShapes(nSeries);
       setupSymbolColors(nSeries);

       // create symbol settings; one per attribute category
       boolean isVisible = true;   // initially all categories are visible
       igvSymbolSettingsList = new ArrayList<SymbolSettings>();

        // plot series uses currently defined series shapes and colors
        for(int series = 0; series < nSeries; ++series) {

            symbolSettings = new SymbolSettings(series, categories[series],
                    igvSymbolShapeList.get(series), igvSymbolColorList.get(series),
                    igvShapeOutlined, igvShapeFilled, isVisible);

            igvSymbolSettingsList.add(symbolSettings);
        }

    }

    /*
    *   Method sets up the shape array for inclusion in symbol settings
    *
    *   Parameters:
    *       nSeries - number of series shapes required.
    * */
    private void setupSymbolShapes(int nSeries) {
        Shape shape;

        // construct the shape array
        igvSymbolShapeList = new ArrayList<Shape>();

         // check if series shapes are the same default
         if(igvDefaultSymbolShape.compareTo(IGVSymbolShape.SeriesSpecific) != 0){

             shape = createSymbolShape(igvDefaultSymbolShape, igvShapeWidth,
                     igvShapeHeight);

             for(int series = 0; series < nSeries; ++series){
                 igvSymbolShapeList.add(shape);
             }
         }
         // series shapes are different - starting with 1st default shape
         // Note: IGVSymbolShape[series index = 0] for IGVSymbolShape.SeriesSpecific
         else {
             IGVSymbolShape[] shapes = IGVSymbolShape.values();


             for(int series = 0; series < nSeries; ++series){
                 int shapesIndex = Math.min(series+1, shapes.length - 1);
                 igvSymbolShapeList.add(createSymbolShape(shapes[shapesIndex], igvShapeWidth, igvShapeHeight));
             }
         }

    }

    /*
    *   Method sets up the color array for inclusion in symbol settings
    *
    *   Parameters:
    *       nSeries - number of series shapes required.
    * */
    private void setupSymbolColors(int nSeries) {

        // construct the color array
        igvSymbolColorList = new ArrayList<Color>();

        for(int series = 0; series < nSeries; ++series){
            int colorIndex = Math.min(series, defaultColor.length - 1);
            igvSymbolColorList.add(defaultColor[colorIndex]);
        }
    }

}
