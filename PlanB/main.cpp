#include <iostream>
#include <queue>
#include <stack>
#include <cassert>
#include <fstream>

#include "gdal_priv.h"
#include "cpl_conv.h"



// A structure that links to a single cell in a Raster
struct RasterCell {
    int x, y; // row and column of the cell
    int elevation;
    int insertion_order;
    bool add_to_list;
    int flow_dir;
    int accumulation = 0;

    // Defines a new link to a cell
    RasterCell(int x, int y, int elevation, int insert_order, int dir) {
        this->x = x;
        this->y = y;
        this->elevation = elevation;
        this->insertion_order = insert_order;
        this->flow_dir = dir;
    }

    // Define the order of the linked cells (to be used in a priority_queue)
    bool operator<(const RasterCell& other) const {
        // to do with statements like if (this->elevation > other.elevation) return false/true;
        if (this->elevation == other.elevation)
        {
            if (this->insertion_order < other.insertion_order)
                return false;
            else if (this->insertion_order > other.insertion_order)
                return true;
        }
        else if (this->elevation <= other.elevation)
            return false;
        else if (this->elevation > other.elevation)
            return true;
       
    }
};

// Storage and access of a raster of a given size
struct Raster {
    std::vector<unsigned int> pixels; // where everything is stored
    std::vector<int> visiting; // store information about if the cell is visited
    std::vector<int> in_queue; // store information about if the cell is added to the priority queue


    int max_x, max_y; // number of columns and rows
    int direction;
    

    // Initialise a raster with x columns and y rows
    Raster(int x, int y) {
        max_x = x;
        max_y = y;
        unsigned int total_pixels = x * y;
        pixels.reserve(total_pixels);
        visiting.reserve(total_pixels);
        
    }

    // Fill values of an entire row
    void add_scanline(const unsigned int* line) {
        for (int i = 0; i < max_x; ++i)
            pixels.push_back(line[i]);
    }

    // Fill entire raster with zeros
    void fill() {
        unsigned int total_pixels = max_x * max_y;
        for (int i = 0; i < total_pixels; ++i)
        {
            pixels.push_back(0);
        }
    }

    // Fill the entire raster with 0 in visiting and in_queue vector
    void fill_visit() {
        unsigned int total_pixels = max_x * max_y;
        for (int i = 0; i < total_pixels; ++i)
        {
            visiting.push_back(0);
            in_queue.push_back(0);
        }
    }


    // Access the value of a raster cell to read or write it
    unsigned int& operator()(int x, int y) {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        return pixels[x + y * max_x];
    }

    // Access the value of a raster cell to read it
    unsigned int operator()(int x, int y) const {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        return pixels[x + y * max_x];
    }

    // Add pixel value for output raster
    void add_value(int x1, int y1, int value) {
        assert(x1 >= 0 && x1 < max_x);
        assert(y1 >= 0 && y1 < max_y);
        pixels[x1 + y1 * max_x] = value;
        
    }

    //change the status from unvisited to visited
    void Is_Visited(int x1, int y1) {
        visiting[x1 + y1 * max_x] = 1;
    }

    // return the status of the pixel whether visited.
    int If_Visited(int x1, int y1) {
        return visiting[x1 + y1 * max_x];
    }

    // update the status after add to the cell_to_process_flow queue
    void Add_to_queue(int x, int y) {
        in_queue[x + y * max_x] = 1;
    }

    // check the status whether add to the cell_to_process_flow queue
    int If_add_to_queue(int x, int y) {
        return in_queue[x + y * max_x];
    }

    void output_accumulation(int& current_line, unsigned int* line)
    {
        for (int i = 0; i < max_y; ++i) 
            line[i] = pixels[i + current_line * max_y];
    }
};

// global variant insert order
int insert_order = 0;

// Write the values in a linked raster cell (useful for debugging)
std::ostream& operator<<(std::ostream& os, const RasterCell& c) {
    os << "{h=" << c.elevation << ", o=" << c.insertion_order << ", x=" << c.x << ", y=" << c.y << "}";
    return os;
};

// write raster result into the tiff file
void output_tiff(std::string filename, Raster input_raster,int max_x, int max_y) {

    int nXSize = max_x;
    int nYSize = max_y;
    double geo_transform[6];
    std::string tiffname = filename;
    GDALAllRegister();
    CPLPushErrorHandler(CPLQuietErrorHandler);
    GDALDataset* gribDataset, * geotiffDataset;
    GDALDriver* driverGeotiff;
    GDALRasterBand* geotiffBand;
    
    gribDataset = (GDALDataset*)GDALOpen("N36E076.hgt", GA_ReadOnly);
    gribDataset->GetGeoTransform(geo_transform);

    driverGeotiff = GetGDALDriverManager()->GetDriverByName("GTiff");
    geotiffDataset = driverGeotiff->Create(tiffname.c_str(), nXSize, nYSize, 1, GDT_Int32, NULL);
    geotiffDataset->SetGeoTransform(geo_transform);
    geotiffDataset->SetProjection(gribDataset->GetProjectionRef());
    unsigned int* rowBuff = (unsigned int*)CPLMalloc(sizeof(unsigned int) * nYSize);

    for (int j = 0; j < nYSize; j++) {
        for (int i = 0; i < nXSize; i++) {
            rowBuff[i] = (unsigned int)input_raster(i, j);
        }
        geotiffDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, j, nYSize, 1, rowBuff, nYSize, 1, GDT_Int32, 0, 0);
    }

    GDALClose(gribDataset);
    GDALClose(geotiffDataset);
    GDALDestroyDriverManager();
}


int main(int argc, const char* argv[]) {
    // Open dataset
    GDALDataset* input_dataset;
    GDALAllRegister();
    input_dataset = (GDALDataset*)GDALOpen("N36E076.hgt", GA_ReadOnly);
    if (input_dataset == NULL) {
        std::cerr << "Couldn't open file" << std::endl;
        return 1;
    }

    // Print dataset info
    double geo_transform[6];
    std::cout << "Driver: " << input_dataset->GetDriver()->GetDescription() << "/" << input_dataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME) << std::endl;;
    std::cout << "Size is " << input_dataset->GetRasterXSize() << "x" << input_dataset->GetRasterYSize() << "x" << input_dataset->GetRasterCount() << std::endl;
    if (input_dataset->GetProjectionRef() != NULL) std::cout << "Projection is '" << input_dataset->GetProjectionRef() << "'" << std::endl;
    if (input_dataset->GetGeoTransform(geo_transform) == CE_None) {
        std::cout << "Origin = (" << geo_transform[0] << ", " << geo_transform[3] << ")" << std::endl;
        std::cout << "Pixel Size = (" << geo_transform[1] << ", " << geo_transform[5] << ")" << std::endl;
    }

    // Print Band 1 info
    GDALRasterBand* input_band;
    int nBlockXSize, nBlockYSize;
    int bGotMin, bGotMax;
    double adfMinMax[2];
    input_band = input_dataset->GetRasterBand(1);
    input_band->GetBlockSize(&nBlockXSize, &nBlockYSize);
    std::cout << "Band 1 Block=" << nBlockXSize << "x" << nBlockYSize << " Type=" << GDALGetDataTypeName(input_band->GetRasterDataType()) << " ColorInterp=" << GDALGetColorInterpretationName(input_band->GetColorInterpretation()) << std::endl;
    adfMinMax[0] = input_band->GetMinimum(&bGotMin);
    adfMinMax[1] = input_band->GetMaximum(&bGotMax);
    if (!(bGotMin && bGotMax))
        GDALComputeRasterMinMax((GDALRasterBandH)input_band, TRUE, adfMinMax);
    std::cout << "Min=" << adfMinMax[0] << " Max=" << adfMinMax[1] << std::endl;

    // Read Band 1 line by line
    int nXSize = input_band->GetXSize();
    int nYSize = input_band->GetYSize();
    Raster input_raster(nXSize, nYSize);
    for (int current_scanline = 0; current_scanline < nYSize; ++current_scanline) {
        unsigned int* scanline = (unsigned int*)CPLMalloc(sizeof(unsigned int) * nXSize);
        if (input_band->RasterIO(GF_Read, 0, current_scanline, nXSize, 1,
            scanline, nXSize, 1, GDT_Int32,
            0, 0) != CPLE_None) {
            std::cerr << "Couldn't read scanline " << current_scanline << std::endl;
            return 1;
        } input_raster.add_scanline(scanline);
        CPLFree(scanline);
    }

    std::cout << "Created raster: " << input_raster.max_x << "x" << input_raster.pixels.size() / input_raster.max_y << " = " << input_raster.pixels.size() << std::endl;
    Raster flow_direction(input_raster.max_x, input_raster.max_y);
    flow_direction.fill();
    flow_direction.fill_visit();

    std::priority_queue<RasterCell, std::deque<RasterCell>> cells_to_process_flow;
    std::deque<RasterCell> cells_to_process_accumulation;

    //add the cell on the boundary(the first and the last row) and choose the initial pixel
    for (int i = 0; i < nXSize; i++)
    {
        int elevation1 = input_raster(i, 0);
        cells_to_process_flow.push(RasterCell(i, 0, elevation1, insert_order, 0));
        insert_order++;
        flow_direction.Add_to_queue(i, 0);

        int elevation2 = input_raster(i, nXSize - 1);
        cells_to_process_flow.push(RasterCell(i, nYSize - 1, elevation2, insert_order, 0));
        flow_direction.Add_to_queue(i, nYSize - 1);
        insert_order++;
    }

    //add the cell on the boundary(the first and the last column) and choose the initial pixel
    for (int j = 1; j < nYSize - 1; j++)
    {
        int elevation3 = input_raster(0, j);
        cells_to_process_flow.push(RasterCell(0, j, elevation3, insert_order, 0));
        flow_direction.Add_to_queue(0, j);
        insert_order++;

        int elevation4 = input_raster(nXSize - 1, j);
        cells_to_process_flow.push(RasterCell(nXSize - 1, j, elevation4, insert_order, 0));
        flow_direction.Add_to_queue(nXSize - 1, j);
        insert_order++;
    }

    //select the lowest elevation cell in the priority queue
    RasterCell start_raster = cells_to_process_flow.top();

    while (cells_to_process_flow.empty() != true)
    {
        // calculate neighbours direction for top_left cell
        start_raster = cells_to_process_flow.top();
        if (start_raster.x == 0 && start_raster.y == 0) {
            if (flow_direction(start_raster.x + 1, start_raster.y) == 0) 
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y) == 0) 
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y, input_raster(start_raster.x + 1, start_raster.y), insert_order, 40));
                    insert_order++;
                } 
                flow_direction(start_raster.x + 1, start_raster.y) = 40;
            }
            if (flow_direction(start_raster.x, start_raster.y + 1) == 0) 
            {
                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y + 1) == 0) {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y + 1, input_raster(start_raster.x, start_raster.y + 1), insert_order, 20));
                    insert_order++;
                } 
                flow_direction(start_raster.x, start_raster.y + 1) = 20;
            }
            if (flow_direction(start_raster.x + 1, start_raster.y + 1) == 0) 
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y + 1) == 0) {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y + 1, input_raster(start_raster.x + 1, start_raster.y + 1), insert_order, 10));
                    insert_order++;
                } 
                flow_direction(start_raster.x + 1, start_raster.y + 1) = 10;
            }
        }

        // calculate neighbours direction for center cells
        else if (start_raster.x > 0 && start_raster.x < nXSize - 1 && start_raster.y>0 && start_raster.y < nYSize - 1)
        {
            if (flow_direction(start_raster.x + 1, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y + 1, input_raster(start_raster.x + 1, start_raster.y + 1), insert_order, 10));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y + 1) = 10;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y + 1, input_raster(start_raster.x - 1, start_raster.y + 1), insert_order, 30));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y + 1) = 30;
            }
            if (flow_direction(start_raster.x, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y + 1, input_raster(start_raster.x, start_raster.y + 1), insert_order, 20));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y + 1) = 20;
            }
            if (flow_direction(start_raster.x + 1, start_raster.y) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y, input_raster(start_raster.x + 1, start_raster.y), insert_order, 40));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y) = 40;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y) == 0) {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y, input_raster(start_raster.x - 1, start_raster.y), insert_order, 50));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y) = 50;
            }
            if (flow_direction(start_raster.x + 1, start_raster.y - 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y - 1, input_raster(start_raster.x + 1, start_raster.y - 1), insert_order, 60));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y - 1) = 60;
            }

            if (flow_direction(start_raster.x, start_raster.y - 1) == 0) {
                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y - 1, input_raster(start_raster.x, start_raster.y - 1), insert_order, 70));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y - 1) = 70;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y - 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y - 1, input_raster(start_raster.x - 1, start_raster.y - 1), insert_order, 80));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y - 1) = 80;
            }

        }
        // calculate neighbours direction for top_right cell
        else if (start_raster.x == nXSize - 1 && start_raster.y == 0)
        {
            if (flow_direction(start_raster.x - 1, start_raster.y) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y, input_raster(start_raster.x - 1, start_raster.y), insert_order, 50));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y) = 50;
            }
            if (flow_direction(start_raster.x, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y + 1, input_raster(start_raster.x, start_raster.y + 1), insert_order, 20));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y + 1) = 20;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y + 1, input_raster(start_raster.x - 1, start_raster.y + 1), insert_order, 30));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y + 1) = 30;
            }

        }
        // calculate neighbours direction for bottom left cell
        else if (start_raster.x == 0 && start_raster.y == nYSize - 1)
        {
            if (flow_direction(start_raster.x, start_raster.y - 1) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y - 1, input_raster(start_raster.x, start_raster.y - 1), insert_order, 70));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y - 1) = 70;
            }

            if (flow_direction(start_raster.x + 1, start_raster.y) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y, input_raster(start_raster.x + 1, start_raster.y), insert_order, 40));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y) = 40;
            }
            if (flow_direction(start_raster.x + 1, start_raster.y - 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y - 1, input_raster(start_raster.x + 1, start_raster.y - 1), insert_order, 60));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y - 1) = 60;
            }

        }
        // calculate neighbours direction for bottom right cell
        else if (start_raster.x == nXSize - 1 && start_raster.y == nYSize - 1)
        {
            if (flow_direction(start_raster.x - 1, start_raster.y) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y, input_raster(start_raster.x - 1, start_raster.y), insert_order, 50));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y) = 50;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y - 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y - 1, input_raster(start_raster.x - 1, start_raster.y - 1), insert_order, 80));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y - 1) = 80;
            }
            if (flow_direction(start_raster.x, start_raster.y - 1) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y - 1, input_raster(start_raster.x, start_raster.y - 1), insert_order, 70));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y - 1) = 70;
            }
        }
        // calculate neighbours direction for left boundary
        else if (start_raster.x == 0 && start_raster.y > 0 && start_raster.y < nYSize - 1)
        {         
            if (flow_direction(start_raster.x + 1, start_raster.y) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y, input_raster(start_raster.x + 1, start_raster.y), insert_order, 40));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y) = 40;
            }
            
            if (flow_direction(start_raster.x, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y + 1, input_raster(start_raster.x, start_raster.y + 1), insert_order, 20));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y + 1) = 20;
            }
            if (flow_direction(start_raster.x + 1, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y + 1, input_raster(start_raster.x + 1, start_raster.y + 1), insert_order, 10));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y + 1) = 10;
            }
            if (flow_direction(start_raster.x, start_raster.y - 1) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y - 1, input_raster(start_raster.x, start_raster.y - 1), insert_order, 70));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y - 1) = 70;
            }
            if (flow_direction(start_raster.x + 1, start_raster.y - 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y - 1, input_raster(start_raster.x + 1, start_raster.y - 1), insert_order, 60));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y - 1) = 60;
            }

        }
        // calculate neighbours direction for right boundary cells
        else if (start_raster.x == nXSize - 1 && start_raster.y > 0 && start_raster.y < nYSize - 1)
        {
            if (flow_direction(start_raster.x, start_raster.y - 1) == 0) {
                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y - 1, input_raster(start_raster.x, start_raster.y - 1), insert_order, 70));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y - 1) = 70;
            }
            if (flow_direction(start_raster.x, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y + 1, input_raster(start_raster.x, start_raster.y + 1), insert_order, 20));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y + 1) = 20;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y, input_raster(start_raster.x - 1, start_raster.y), insert_order, 50));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y) = 50;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y - 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y - 1, input_raster(start_raster.x - 1, start_raster.y - 1), insert_order, 80));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y - 1) = 80;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y + 1, input_raster(start_raster.x - 1, start_raster.y + 1), insert_order, 30));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y + 1) = 30;
            }

        }
        // calculate neighbours direction top boundary cells
        else if (start_raster.x > 0 && start_raster.x < nXSize - 1 && start_raster.y == 0)
        {
            if (flow_direction(start_raster.x + 1, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y + 1, input_raster(start_raster.x + 1, start_raster.y + 1), insert_order, 10));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y + 1) = 10;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y + 1, input_raster(start_raster.x - 1, start_raster.y + 1), insert_order, 30));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y + 1) = 30;
            }
            if (flow_direction(start_raster.x, start_raster.y + 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y + 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y + 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y + 1, input_raster(start_raster.x, start_raster.y + 1), insert_order, 20));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y + 1) = 20;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y, input_raster(start_raster.x - 1, start_raster.y), insert_order, 50));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y) = 50;
            }
            if (flow_direction(start_raster.x + 1, start_raster.y) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y, input_raster(start_raster.x + 1, start_raster.y), insert_order, 40));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y) = 40;
            }
        }
        // calculate neighbours direction for bottom boundary cell
        else if (start_raster.x > 0 && start_raster.x < nXSize - 1 && start_raster.y == nYSize - 1)
        {
            if (flow_direction(start_raster.x + 1, start_raster.y) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y, input_raster(start_raster.x + 1, start_raster.y), insert_order, 40));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y) = 40;
            }
            if (flow_direction(start_raster.x - 1, start_raster.y) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y, input_raster(start_raster.x - 1, start_raster.y), insert_order, 50));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y) = 50;
            }
            if (flow_direction(start_raster.x + 1, start_raster.y - 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x + 1, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x + 1, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x + 1, start_raster.y - 1, input_raster(start_raster.x + 1, start_raster.y - 1), insert_order, 60));
                    insert_order++;
                }
                flow_direction(start_raster.x + 1, start_raster.y - 1) = 60;
            }

            if (flow_direction(start_raster.x, start_raster.y - 1) == 0) {

                if (flow_direction.If_add_to_queue(start_raster.x, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x, start_raster.y - 1, input_raster(start_raster.x, start_raster.y - 1), insert_order, 70));
                    insert_order++;
                }
                flow_direction(start_raster.x, start_raster.y - 1) = 70;
            }

            if (flow_direction(start_raster.x - 1, start_raster.y - 1) == 0)
            {
                if (flow_direction.If_add_to_queue(start_raster.x - 1, start_raster.y - 1) == 0)
                {
                    flow_direction.Add_to_queue(start_raster.x - 1, start_raster.y - 1);
                    cells_to_process_flow.push(RasterCell(start_raster.x - 1, start_raster.y - 1, input_raster(start_raster.x - 1, start_raster.y - 1), insert_order, 80));
                    insert_order++;
                }
                flow_direction(start_raster.x - 1, start_raster.y - 1) = 80;
            }

        }

        // change the visit status of each cell
        flow_direction.Is_Visited(start_raster.x, start_raster.y);
        // add to the cell to the stack to calculate flow accumulation later
        cells_to_process_accumulation.push_back(start_raster);
        //delete and select new center cell
        cells_to_process_flow.pop();
    }
    // output tif file
    output_tiff("flow_direction.tif", flow_direction, nXSize, nYSize);
    std::cout << "finish output flow direction tiff file" << std::endl;

    // initialize the flow accumulation raster
    Raster flow_accumulation(nXSize, nYSize);
    flow_accumulation.fill();
    RasterCell cell = cells_to_process_accumulation.back();
    //iterate and calculate flow accumulation base on the flow direction
    while (cells_to_process_accumulation.empty() != true)
    {
        cell = cells_to_process_accumulation.back();
        int x = cell.x;
        int y = cell.y;
        int dir = flow_direction(x,y);
        switch (dir)
        {
        case 10:
            if (x - 1 >= 0 && x - 1 < nXSize && y - 1 >= 0 && y - 1 < nYSize) 
            {
                flow_accumulation(x - 1, y - 1) += flow_accumulation(x, y) + 1;
                break;
            }
        case 20:
            if (x >= 0 && x < nXSize && y - 1 >= 0 && y - 1 < nYSize) 
            {
                flow_accumulation(x, y - 1) += flow_accumulation(x, y) + 1;
                break;
            }
        case 30:
            if (x + 1 >= 0 && x + 1 < nXSize && y - 1 >= 0 && y - 1 < nYSize) 
            {
                flow_accumulation(x + 1, y - 1) += flow_accumulation(x, y) + 1;
                break;
            }
        case 40:
            if (x - 1 >= 0 && x - 1 < nXSize && y >= 0 && y < nYSize) 
            {
                flow_accumulation(x - 1, y) += flow_accumulation(x, y) + 1;
                break;
            }
        case 50:
            if (x + 1 >= 0 && x + 1 < nXSize && y >= 0 && y < nYSize) 
            {
                flow_accumulation(x + 1, y) += flow_accumulation(x, y) + 1;
                break;
            }
        case 60:
            if (x - 1 >= 0 && x - 1 < nXSize && y + 1 >= 0 && y + 1 < nYSize) {
                flow_accumulation(x - 1, y + 1) += flow_accumulation(x, y) + 1;
                break;
            }
        case 70:
            if (x >= 0 && x < nXSize && y + 1 >= 0 && y + 1 < nYSize) {
                flow_accumulation(x, y + 1) += flow_accumulation(x, y) + 1;
                break;
            }
        case 80:
            if (x + 1 >= 0 && x + 1 < nXSize && y + 1 >= 0 && y + 1 < nYSize) {
                flow_accumulation(x + 1, y + 1) += flow_accumulation(x, y) + 1;
                break;
            }
        case 0:
            if (x >= 0 && x < nXSize && y >= 0 && y < nYSize) {
                flow_accumulation(x, y) = 0;
                break;
            }
        }
        cells_to_process_accumulation.pop_back();
    }
    //output flow accumulation raster
    output_tiff("flow_accumulation.tif", flow_accumulation, nXSize, nYSize);
    std::cout << "finish output flow_accumulation tiff file" << std::endl;
    GDALClose(input_dataset);
    return 0;
}
