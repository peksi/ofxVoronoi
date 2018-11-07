#pragma once

// openFrameworks
#include "ofMain.h"

class ofxVoronoiCell {
  public:
    vector<glm::vec2> pts;
    glm::vec2 pt;
};

class ofxVoronoi {
private:
    ofRectangle bounds;
    vector<glm::vec2> points;
    vector<ofxVoronoiCell> cells;
    
public:
    ofxVoronoi();
    ~ofxVoronoi();
    
    void clear();
    void generate(bool ordered=true);
    vector<ofPolyline> draw(bool curve);
    
    bool isBorder(glm::vec2 _pt);
    
    void setBounds(ofRectangle _bounds);
    void setPoints(vector<glm::vec2> _points);
    void addPoint(glm::vec2 _point);
    void addPoints(vector<glm::vec2> _points);
    
    ofRectangle getBounds();
    vector<glm::vec2>& getPoints();
    vector <ofxVoronoiCell>& getCells();
    ofxVoronoiCell& getCell(glm::vec2 _point, bool approximate=false);
    
    //borg
    void relax();
};
