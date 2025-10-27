#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <string>

class Image{
private: 
    int width;
    int height;
    std::vector<unsigned int> pixels;
    
public:
    Image(int w, int h);    // constructor for write
    Image(const std::string& filename);     // constructor for read 

    int getWidth() const { return width; }
    int getHeight() const { return height; }
    
    int clamp(int value, int min, int max);

    // Pixel operations
    void setPixel(int x, int y, int r, int g, int b);
    void getPixel(int x, int y, int& r, int& g, int& b) const;

    // File I/O
    void write(const std::string& filename) const;
    void read(const std::string& filename);

private:
    int getIndex(int x, int y) const;

};

#endif