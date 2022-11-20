#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <Utils.h>
#include <algorithm>    // sort()
#include <ctime>        // time()
#include <sstream>      // stringstream
#include <cmath>        // pow()
#include <map>

#define WIDTH 360
#define HEIGHT 360
#define epsilon 0.00001

uint32_t backgroundColour = (255 << 24) + (255 << 16) + (255 << 8) + 255;
glm::vec3 cam_pos = glm::vec3(0.0f, 0.0f, 12.0f);
glm::mat3 cam_orientation = glm::mat3{
	glm::vec3(1.0f, 0.0f, 0.0f),
	glm::vec3(0.0f, 1.0f, 0.0f),
	glm::vec3(0.0f, 0.0f, 1.0f),
};
glm::vec3 light = glm::vec3(1.0f, 1.0f, 3.0f);
//glm::vec3 light = glm::vec3(0.5f, 0.5f, 3.0f);
float lightIntensity = 2.0f;
float brightnessThreshold = 0.2f;
float focal = 2.0f;

std::vector<std::vector<float>> depthMap(HEIGHT, std::vector<float>(WIDTH, 0.0f));

struct Material {
	Colour c;
	bool texturePresent = false;
	TextureMap source;
};

struct LogicallyNamedPoints {
	CanvasPoint top;
	CanvasPoint left;
	CanvasPoint right;
	CanvasPoint bottom;
};

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	std::vector<float> v;
	if ((numberOfValues == 1) || (numberOfValues == 0)) return v;
	for (int i = 0; i < numberOfValues; i++)
		v.push_back(from + i * (to - from) / (numberOfValues - 1));
	return v;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> v;
	if ((numberOfValues == 1) || (numberOfValues == 0)) return v;
	for (int i = 0; i < numberOfValues; i++) {
		float k = i;
		float n = numberOfValues;
		v.push_back(from + k * (to - from) / (n - 1));
	}
	return v;
}

float depth(float l_x, float r_x, float p_x, float l_d, float r_d) {
	if (l_d == 0.0f || r_d == 0.0f) return 0.0f;
	return 1 / (((r_x - p_x) / (r_x - l_x)) * ((1 / l_d) - (1 / r_d)) + (1 / r_d));
}

void drawLine(DrawingWindow& window, CanvasPoint from, CanvasPoint to, Colour colour) {
	uint32_t setColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;

	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;

	float noOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff / noOfSteps;
	float yStepSize = yDiff / noOfSteps;

	for (float i = 0.0; i < noOfSteps; i++) {
		float x = from.x + (xStepSize * i);
		float y = from.y + (yStepSize * i);
		float d = depth(from.x, to.x, x, from.depth, to.depth);
		if (x < WIDTH && y < HEIGHT && x >= 0 && y >= 0) {
			if (d >= depthMap[y][x]) {
				window.setPixelColour(x, y, setColour);
				depthMap[y][x] = d;
			}
		}
	}
}

void drawStrokedTriangle(DrawingWindow& window, CanvasTriangle t, Colour colour) {
	drawLine(window, t.vertices[0], t.vertices[1], colour);
	drawLine(window, t.vertices[1], t.vertices[2], colour);
	drawLine(window, t.vertices[2], t.vertices[0], colour);
}

bool compare(CanvasPoint i, CanvasPoint j) { return (i.y < j.y); }

void drawFilledTriangle(DrawingWindow& window, CanvasTriangle t, ModelTriangle mt) {
	drawStrokedTriangle(window, t, mt.colour);
	uint32_t setColour = (255 << 24) + (mt.colour.red << 16) + (mt.colour.green << 8) + mt.colour.blue;

	sort(t.vertices.begin(), t.vertices.end(), compare);

	LogicallyNamedPoints points;
	points.top = t.vertices[0];
	points.bottom = t.vertices[2];
	CanvasPoint middle = t.vertices[1];

	float ic_x = points.top.x - (points.top.y - middle.y) / ((points.top.y - points.bottom.y) / (points.top.x - points.bottom.x));
	CanvasPoint ic = CanvasPoint(ic_x, middle.y, depth(points.top.x, points.bottom.x, ic_x, points.top.depth, points.bottom.depth));

	if (ic.x > middle.x) {
		points.right = ic;
		points.left = middle;
	}
	else {
		points.right = middle;
		points.left = ic;
	}

	std::vector<float> topLeft = interpolateSingleFloats(points.top.x, points.left.x, abs(points.top.y - middle.y) + 1);
	std::vector<float> topRight = interpolateSingleFloats(points.top.x, points.right.x, abs(points.top.y - middle.y) + 1);
	std::vector<float> bottomLeft = interpolateSingleFloats(points.left.x, points.bottom.x + 1, abs(middle.y - points.bottom.y) + 2);
	std::vector<float> bottomRight = interpolateSingleFloats(points.right.x, points.bottom.x + 1, abs(middle.y - points.bottom.y) + 2);

	std::vector<float> tl_depth;
	std::vector<float> tr_depth;
	for (int i = 0; i < topLeft.size(); i++) {
		tl_depth.push_back(depth(points.top.x, points.left.x, topLeft[i], points.top.depth, points.left.depth));
		tr_depth.push_back(depth(points.top.x, points.right.x, topRight[i], points.top.depth, points.right.depth));
	}
	std::vector<float> bl_depth;
	std::vector<float> br_depth;
	for (int i = 0; i < bottomLeft.size(); i++) {
		bl_depth.push_back(depth(points.left.x, points.bottom.x + 1, bottomLeft[i], points.left.depth, points.bottom.depth));
		br_depth.push_back(depth(points.right.x, points.bottom.x + 1, bottomRight[i], points.right.depth, points.bottom.depth));
	}

	int i = 0; int k = 0;
	for (float y = points.top.y; i < topLeft.size(); y++, i++) {
		for (float x = topLeft[i]; x < topRight[i]; x++) {
			if (x < WIDTH && y < HEIGHT && x >= 0 && y >= 0) {
				float d = depth(topLeft[i], topRight[i], x, tl_depth[i], tr_depth[i]);
				if (d >= depthMap[y][x]) {
					window.setPixelColour(x, y, setColour);
					depthMap[y][x] = d;
				}
			}
		}
	}
	for (float y = middle.y; k < bottomLeft.size(); y++, k++) {
		for (float x = bottomLeft[k]; x < bottomRight[k]; x++) {
			if (x < WIDTH && y < HEIGHT && x >= 0 && y >= 0) {
				float d = depth(bottomLeft[k], bottomRight[k], x, bl_depth[k], br_depth[k]);
				if (d >= depthMap[y][x]) {
					window.setPixelColour(x, y, setColour);
					depthMap[y][x] = d;
				}
			}
		}
	}
}

uint32_t getPerspectivePixel(TextureMap source, int n, float x_l, float x_r, int length, float z0, float z1, float c0, float c1, float q) {
	float c = (c0 / z0 * (1 - q) + c1 / z1 * q) / (1 / z0 * (1 - q) + (1 / z1 * q));
	std::vector<float> xs = interpolateSingleFloats(x_l, x_r, length + 1);
	float x_s = xs[n];
	uint32_t colour = source.pixels[int(c) * source.width + int(x_s)];

	return colour;
}

void drawTextureTriangle(DrawingWindow& window, TextureMap source, CanvasTriangle t) {
	sort(t.vertices.begin(), t.vertices.end(), compare);

	LogicallyNamedPoints points;
	points.top = t.vertices[0];
	points.bottom = t.vertices[2];
	CanvasPoint middle = t.vertices[1];

	float ic_x = points.top.x - (points.top.y - middle.y) / ((points.top.y - points.bottom.y) / (points.top.x - points.bottom.x));
	CanvasPoint ic = CanvasPoint(ic_x, middle.y, depth(points.top.x, points.bottom.x, ic_x, points.top.depth, points.bottom.depth));
	float ratio = (ic.y - points.top.y) / (points.bottom.y - points.top.y);
	ic.texturePoint.x = points.top.texturePoint.x + (points.bottom.texturePoint.x - points.top.texturePoint.x) * ratio;
	ic.texturePoint.y = points.top.texturePoint.y + (points.bottom.texturePoint.y - points.top.texturePoint.y) * ratio;
	ic.texturePoint = TexturePoint(ic.texturePoint.x, ic.texturePoint.y);

	if (ic.x > middle.x) {
		points.right = ic;
		points.left = middle;
	}
	else {
		points.right = middle;
		points.left = ic;
	}

	std::vector<float> topLeft = interpolateSingleFloats(points.top.x, points.left.x, abs(points.top.y - middle.y) + 1);
	std::vector<float> topRight = interpolateSingleFloats(points.top.x, points.right.x, abs(points.top.y - middle.y) + 1);
	std::vector<float> bottomLeft = interpolateSingleFloats(points.left.x, points.bottom.x + 1, abs(middle.y - points.bottom.y) + 2);
	std::vector<float> bottomRight = interpolateSingleFloats(points.right.x, points.bottom.x + 1, abs(middle.y - points.bottom.y) + 2);
	// interpolate the source
	std::vector<float> topLeft_s = interpolateSingleFloats(points.top.texturePoint.x, points.left.texturePoint.x, abs(points.top.y - middle.y) + 1);
	std::vector<float> topRight_s = interpolateSingleFloats(points.top.texturePoint.x, points.right.texturePoint.x, abs(points.top.y - middle.y) + 1);
	std::vector<float> source_y_tl = interpolateSingleFloats(points.top.texturePoint.y, points.left.texturePoint.y, abs(middle.y - points.top.y) + 1);
	std::vector<float> source_y_tr = interpolateSingleFloats(points.top.texturePoint.y, points.right.texturePoint.y, abs(middle.y - points.top.y) + 1);
	std::vector<float> bottomLeft_s = interpolateSingleFloats(points.left.texturePoint.x, points.bottom.texturePoint.x, abs(middle.y - points.bottom.y) + 2);
	std::vector<float> bottomRight_s = interpolateSingleFloats(points.right.texturePoint.x, points.bottom.texturePoint.x, abs(middle.y - points.bottom.y) + 2);
	std::vector<float> source_y_bl = interpolateSingleFloats(points.left.texturePoint.y, points.bottom.texturePoint.y, abs(middle.y - points.bottom.y) + 2);
	std::vector<float> source_y_br = interpolateSingleFloats(points.right.texturePoint.y, points.bottom.texturePoint.y, abs(middle.y - points.bottom.y) + 2);

	std::vector<float> tl_depth;
	std::vector<float> tr_depth;
	for (int i = 0; i < topLeft.size(); i++) {
		tl_depth.push_back(depth(points.top.x, points.left.x, topLeft[i], points.top.depth, points.left.depth));
		tr_depth.push_back(depth(points.top.x, points.right.x, topRight[i], points.top.depth, points.right.depth));
	}
	std::vector<float> bl_depth;
	std::vector<float> br_depth;
	for (int i = 0; i < bottomLeft.size(); i++) {
		bl_depth.push_back(depth(points.left.x, points.bottom.x + 1, bottomLeft[i], points.left.depth, points.bottom.depth));
		br_depth.push_back(depth(points.right.x, points.bottom.x + 1, bottomRight[i], points.right.depth, points.bottom.depth));
	}

	int i = 0; int k = 0; // ith and kth row
	for (float y = points.top.y; i < topLeft.size() && topLeft.size() > 2; y++, i++) {
		int n = 0; // nth point in current row
		for (float x = topLeft[i]; x < topRight[i]; x++, n++) {
			if (x < WIDTH && y < HEIGHT && x >= 0 && y >= 0) {
				float d = depth(topLeft[i], topRight[i], x, tl_depth[i], tr_depth[i]);
				if (d >= depthMap[y][x]) {
					int length = topRight[i] - topLeft[i] + 1;
					uint32_t colour = getPerspectivePixel(source, n, topLeft_s[i], topRight_s[i], length, 1 / points.top.depth,
						1 / middle.depth, points.top.texturePoint.y, middle.texturePoint.y, i / (middle.y - points.top.y));
					window.setPixelColour(x, y, colour);
					depthMap[y][x] = d;
				}
			}
		}
	}
	for (float y = middle.y; k < bottomLeft.size() && bottomLeft.size() > 2; y++, k++) {
		int n = 0;
		for (float x = bottomLeft[k]; x < bottomRight[k]; x++, n++) {
			if (x < WIDTH && y < HEIGHT && x >= 0 && y >= 0) {
				float d = depth(bottomLeft[k], bottomRight[k], x, bl_depth[k], br_depth[k]);
				if (d >= depthMap[y][x]) {
					int length = bottomRight[k] - bottomLeft[k] + 1;
					uint32_t colour = getPerspectivePixel(source, n, bottomLeft_s[k], bottomRight_s[k], length, 1 / middle.depth,
						1 / points.bottom.depth, middle.texturePoint.y, points.bottom.texturePoint.y, k / (points.bottom.y - middle.y));
					window.setPixelColour(x, y, colour);
					depthMap[y][x] = d;
				}
			}
		}
	}
}

void draw(DrawingWindow& window) {
	window.clearPixels();

	glm::vec3 topLeft(0, 0, 0);
	glm::vec3 topRight(0, 0, 0);
	glm::vec3 bottomRight(0, 0, 0);
	glm::vec3 bottomLeft(0, 0, 0);

	std::vector<glm::vec3> tl2bl = interpolateThreeElementValues(topLeft, bottomLeft, HEIGHT + 1);
	std::vector<glm::vec3> tr2br = interpolateThreeElementValues(topRight, bottomRight, HEIGHT + 1);

	for (size_t y = 0; y < window.height; y++) {
		std::vector<glm::vec3> v = interpolateThreeElementValues(tl2bl[y], tr2br[y], WIDTH + 1);
		for (size_t x = 0; x < window.width; x++) {
			uint32_t colour = (255 << 24) + (int(v[x][0]) << 16) + (int(v[x][1]) << 8) + int(v[x][2]);
			window.setPixelColour(x, y, colour);
		}
	}
}

glm::mat4 translate(glm::vec3 v) {
	return {
		glm::vec4(1.0f, 0.0f, 0.0f, 0.0f),
		glm::vec4(0.0f, 1.0f, 0.0f, 0.0f),
		glm::vec4(0.0f, 0.0f, 1.0f, 0.0f),
		glm::vec4(v, 1.0f),
	};
}

glm::mat3 rotateX(float theta) {
	return {
		glm::vec3(1.0f, 0.0f, 0.0f),
		glm::vec3(0.0f, cos(theta), sin(theta)),
		glm::vec3(0.0f, -sin(theta), cos(theta)),
	};
}

glm::mat3 rotateY(float theta) {
	return {
		glm::vec3(cos(theta), 0.0f, -sin(theta)),
		glm::vec3(0.0f, 1.0f, 0.0f),
		glm::vec3(sin(theta), 0.0f, cos(theta)),
	};
}

glm::mat3 rotateZ(float theta) {
	return {
		glm::vec3(cos(theta), sin(theta), 0.0f),
		glm::vec3(-sin(theta), cos(theta), 0.0f),
		glm::vec3(0.0f, 0.0f, 1.0f),
	};
}

void translateCam(glm::vec3 v) {
	cam_pos = glm::vec3(translate(v) * glm::vec4(cam_pos, 1));
}

void rotateCam(glm::vec3 r) {
	cam_pos = rotateX(r.x) * rotateY(r.y) * rotateZ(r.z) * cam_pos;
}

void rotateSource(glm::mat3 m) {
	cam_orientation = cam_orientation * m;
	light = light * cam_orientation;
}

void lookAt() {
	glm::vec3 z = glm::normalize(cam_pos);
	glm::vec3 x = glm::cross({ 0.0f, 1.0f, 0.0f }, z);
	glm::vec3 y = glm::cross(z, x);

	cam_orientation = { x, y, z };
}

void update(DrawingWindow& window) {
	// Function for performing animation (shifting artifacts or moving the camera)
	//srand(time(NULL));
}

void handleEvent(SDL_Event event, DrawingWindow& window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) { std::cout << "LEFT" << std::endl; lookAt(); }
		else if (event.key.keysym.sym == SDLK_RIGHT) { std::cout << "RIGHT" << std::endl; rotateSource(rotateX(-0.2)); }
		else if (event.key.keysym.sym == SDLK_UP) { std::cout << "UP" << std::endl; rotateCam(glm::vec3(0.08, 0.08, 0.08)); }
		else if (event.key.keysym.sym == SDLK_DOWN) { std::cout << "DOWN" << std::endl; light += glm::vec3(0.0, -0.2, 0.0); }
		else if (event.key.keysym.sym == SDLK_F1) { std::cout << "U" << std::endl; translateCam(glm::vec3(0.0, -0.2, 0.0)); }
		else if (event.key.keysym.sym == SDLK_F2) { std::cout << "D" << std::endl; translateCam(glm::vec3(0.0, 0.2, 0.0)); }
		else if (event.key.keysym.sym == SDLK_F3) { std::cout << "F" << std::endl; translateCam(glm::vec3(0.0, 0.0, -0.2)); }
		else if (event.key.keysym.sym == SDLK_F4) { std::cout << "B" << std::endl; translateCam(glm::vec3(0.0, 0.0, 0.2)); }
		else if (event.key.keysym.sym == SDLK_F5) { std::cout << "L" << std::endl; translateCam(glm::vec3(0.2, 0.0, 0.0)); }
		else if (event.key.keysym.sym == SDLK_F6) { std::cout << "R" << std::endl; translateCam(glm::vec3(-0.2, 0.0, 0.0)); }
		else if (event.key.keysym.sym == SDLK_F7) { std::cout << "S" << std::endl; rotateCam(glm::vec3(0.1, 0.1, 0.1)); }
		else if (event.key.keysym.sym == SDLK_F8) { std::cout << "W" << std::endl; rotateCam(glm::vec3(-0.1, -0.1, -0.1)); }
		else if (event.key.keysym.sym == SDLK_F9) { std::cout << "X" << std::endl; rotateSource(rotateX(-0.2)); }
		else if (event.key.keysym.sym == SDLK_F10) { std::cout << "Y" << std::endl; rotateSource(rotateY(0.4)); }
		else if (event.key.keysym.sym == SDLK_F11) { std::cout << "Z" << std::endl; rotateSource(rotateZ(0.2)); }
		else if (event.key.keysym.sym == SDLK_F12) { std::cout << "A" << std::endl; lookAt(); }
	}
	else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

std::map<std::string, Material> loadMtl(const char* filename) {
	std::map<std::string, Material> palette;

	std::ifstream in_file(filename);
	std::string nextline;
	glm::vec3 c_temp;
	std::string n_temp;

	if (!in_file.is_open()) throw "ERROR::OBJLOADER::COULD NOT OPEN.";

	while (std::getline(in_file, nextline)) {
		if (nextline.size() > 0) {
			std::vector<std::string> ws = split(nextline, ' ');
			if (ws[0] == "newmtl") {
				n_temp = ws[1];
				std::getline(in_file, nextline);
				std::vector<std::string> ws = split(nextline, ' ');
				c_temp.x = std::stof(ws[1]) * 255;
				c_temp.y = std::stof(ws[2]) * 255;
				c_temp.z = std::stof(ws[3]) * 255;
				palette[n_temp].c = Colour(n_temp, c_temp.x, c_temp.y, c_temp.z);
			}
			else if (ws[0] == "map_Kd") {
				TextureMap source = TextureMap(ws[1]);
				palette[n_temp].texturePresent = true;
				palette[n_temp].source = source;
			}
		}
	}

	return palette;
}

std::vector<ModelTriangle> loadOBJ(const char* filename, std::map<std::string, Material> palette, float scale) {
	std::ifstream in_file(filename);
	std::string nextline;
	std::stringstream ss;
	std::string prefix;

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> v_texture;
	std::vector<glm::vec3> v_normal;
	std::vector<glm::vec3> facets;
	std::vector<glm::vec3> facets_vt;
	std::vector<glm::vec3> facets_vn;
	std::vector<Colour> colours;
	std::vector<int> c_indices;
	std::vector<ModelTriangle> mts;

	Colour c_temp;
	glm::vec3 v_temp;
	glm::vec2 vt_temp;
	glm::vec3 vn_temp;
	glm::vec3 v_indices;
	glm::vec3 vt_indices;
	glm::vec3 vn_indices;
	int c_indice = 0;

	if (!in_file.is_open()) throw "ERROR::OBJLOADER::COULD NOT OPEN.";

	while (std::getline(in_file, nextline)) {
		if (nextline.size() > 0) {
			ss.clear();
			ss.str(nextline);
			ss >> prefix;

			if (prefix == "v") {
				ss >> v_temp.x >> v_temp.y >> v_temp.z;
				vertices.push_back(v_temp * scale);
			}
			else if (prefix == "vt") {
				ss >> vt_temp.x >> vt_temp.y;
				v_texture.push_back(vt_temp);
			}
			else if (prefix == "vn") {
				ss >> vn_temp.x >> vn_temp.y >> vn_temp.z;
				v_normal.push_back(vn_temp);
			}
			else if (prefix == "f") {
				std::vector<std::string> fs = split(nextline, ' ');

				std::string component = fs[1].substr(fs[1].find_first_of("/"), fs[1].size());
				std::string componentN = fs[1].substr(fs[1].find_last_of("/"), fs[1].size());
				if (component == componentN) {
					vn_indices.x = vn_indices.y = vn_indices.z = -1;
					facets_vn.push_back(vn_indices);
					if (component != "/") {
						vt_indices.x = std::stoi(component.substr(1));
						component = fs[2].substr(fs[2].find_first_of("/"), fs[2].size());
						vt_indices.y = std::stoi(component.substr(1));
						component = fs[3].substr(fs[3].find_first_of("/"), fs[3].size());
						vt_indices.z = std::stoi(component.substr(1));
						facets_vt.push_back(vt_indices);
					}
					else {
						vt_indices.x = vt_indices.y = vt_indices.z = -1;
						facets_vt.push_back(vt_indices);
					}
				}
				else {
					vn_indices.x = std::stoi(componentN.substr(1));
					componentN = fs[2].substr(fs[2].find_last_of("/"), fs[2].size());
					vn_indices.y = std::stoi(componentN.substr(1));
					componentN = fs[3].substr(fs[3].find_last_of("/"), fs[3].size());
					vn_indices.z = std::stoi(componentN.substr(1));
					facets_vn.push_back(vn_indices);

					vt_indices.x = vt_indices.y = vt_indices.z = -1;
					facets_vt.push_back(vt_indices);
				}

				v_indices.x = std::stoi(fs[1]);
				v_indices.y = std::stoi(fs[2]);
				v_indices.z = std::stoi(fs[3]);

				facets.push_back(v_indices);
				c_indice++;
			}
			else if (prefix == "usemtl") {
				std::vector<std::string> ms = split(nextline, ' ');
				c_temp = palette[ms[1]].c;
				c_temp.name = ms[1];
				colours.push_back(c_temp);
				c_indices.push_back(c_indice);
			}
		}
	}

	int i = 0;
	c_indices.push_back(c_indice);
	for (int k = 0; k < facets.size(); k++) {
		if (c_indices.size() > 1 && k == c_indices[1]) {
			c_indices.erase(c_indices.begin());
			i++;
		}
		else colours.push_back(Colour(255, 255, 255));

		ModelTriangle modelTriangle =
			ModelTriangle(vertices[int(facets[k].x) - 1], vertices[int(facets[k].y) - 1], vertices[int(facets[k].z) - 1], colours[i]);

		glm::vec3 e1 = modelTriangle.vertices[1] - modelTriangle.vertices[0];
		glm::vec3 e2 = modelTriangle.vertices[2] - modelTriangle.vertices[1];
		modelTriangle.normal = glm::normalize(glm::cross(e1, e2));

		if (facets_vt[k].x != -1) {
			TexturePoint t1 = TexturePoint(v_texture[int(facets_vt[k].x) - 1].x * palette[colours[i].name].source.width,
				(1 - v_texture[int(facets_vt[k].x) - 1].y) * palette[colours[i].name].source.height);
			TexturePoint t2 = TexturePoint(v_texture[int(facets_vt[k].y) - 1].x * palette[colours[i].name].source.width,
				(1 - v_texture[int(facets_vt[k].y) - 1].y) * palette[colours[i].name].source.height);
			TexturePoint t3 = TexturePoint(v_texture[int(facets_vt[k].z) - 1].x * palette[colours[i].name].source.width,
				(1 - v_texture[int(facets_vt[k].z) - 1].y) * palette[colours[i].name].source.height);
			modelTriangle.texturePoints = { t1, t2, t3 };
		}

		if (facets_vn[k].x != -1) {
			glm::vec3 vn1 = glm::vec3(v_normal[int(facets_vn[k].x) - 1].x, v_normal[int(facets_vn[k].x) - 1].y, v_normal[int(facets_vn[k].x) - 1].z);
			glm::vec3 vn2 = glm::vec3(v_normal[int(facets_vn[k].y) - 1].x, v_normal[int(facets_vn[k].y) - 1].y, v_normal[int(facets_vn[k].y) - 1].z);
			glm::vec3 vn3 = glm::vec3(v_normal[int(facets_vn[k].z) - 1].x, v_normal[int(facets_vn[k].z) - 1].y, v_normal[int(facets_vn[k].z) - 1].z);
			modelTriangle.vn = { vn1, vn2, vn3 };
		}
		else {
			glm::vec3 vn = modelTriangle.normal;
			modelTriangle.vn = { vn, vn, vn };
		}

		mts.push_back(modelTriangle);
	}

	return mts;
}

RayTriangleIntersection getClosestIntersection(glm::vec3 start, glm::vec3 rayDirection, std::vector<ModelTriangle> mts) {
	glm::vec3 solution = glm::vec3(INT_MAX, INT_MAX, INT_MAX);
	int k = 0;

	for (int i = 0; i < mts.size(); i++) {
		glm::vec3 e0 = (mts[i].vertices[1] - mts[i].vertices[0]) * cam_orientation;
		glm::vec3 e1 = (mts[i].vertices[2] - mts[i].vertices[0]) * cam_orientation;
		glm::vec3 SPVector = start - mts[i].vertices[0] * cam_orientation;
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		if (possibleSolution.x < solution.x + epsilon && possibleSolution.x > epsilon) {
			if (possibleSolution.y + epsilon >= 0.0f && possibleSolution.y <= 1.0f + epsilon && possibleSolution.z + epsilon >= 0.0f &&
				possibleSolution.z <= 1.0f + epsilon && possibleSolution.y + possibleSolution.z <= 1.0f + epsilon) {
				solution = possibleSolution;
				k = i;
			}
		}
	}

	RayTriangleIntersection intersection = RayTriangleIntersection(solution, solution.x, mts[k], k);

	return intersection;
}

float incidenceLighting(glm::vec3 vn, glm::vec3 lightD) {
	glm::vec3 normal = vn * cam_orientation;
	float incidence = glm::dot(-lightD, normal);
	return std::max(0.0f, incidence);
}

float specularLighting(glm::vec3 vn, glm::vec3 lightD, glm::vec3 rayD) {
	glm::vec3 normal = vn * cam_orientation;
	glm::vec3 reflection = glm::normalize(lightD - 2.0f * normal * (glm::dot(lightD, normal)));
	float poweredSpecular = glm::dot(reflection, -rayD);
	if (poweredSpecular < 0) return 0.0f;
	float specular = pow(poweredSpecular, 256);
	return specular;
}

float normalShading(std::vector<ModelTriangle> mts, int k, glm::vec3 rayD, float camDistance) {
	glm::vec3 p = rayD * camDistance + cam_pos;
	glm::vec3 lightDirection = glm::normalize(p - light);
	RayTriangleIntersection lightRay = getClosestIntersection(light, lightDirection, mts);
	float brightness;

	if (k == lightRay.triangleIndex) {
		float incidence = incidenceLighting(mts[k].normal, lightDirection);
		float specular = specularLighting(mts[k].normal, lightDirection, rayD);
		float proximity = lightIntensity / lightRay.distanceFromCamera;
		brightness = std::max(specular + incidence * proximity, brightnessThreshold);
	}
	else {
		brightness = brightnessThreshold * 0.9;
	}

	return brightness;
}

float gouraudShading(std::vector<ModelTriangle> mts, int k, glm::vec3 rayD, float camDistance) {
	glm::vec3 p = rayD * camDistance + cam_pos;
	glm::vec3 pLightDirection = glm::normalize(p - light);
	RayTriangleIntersection pLightRay = getClosestIntersection(light, pLightDirection, mts);
	std::vector<float> vBrightness;
	float pBrightness;

	if (k == pLightRay.triangleIndex || mts[k].vn[0] != mts[k].vn[1]) {
		for (int i = 0; i < mts[k].vn.size(); i++) {
			glm::vec3 lightDirection = glm::normalize(mts[k].vertices[i] - light);

			float incidence = incidenceLighting(mts[k].vn[i], lightDirection);
			float specular = specularLighting(mts[k].vn[i], lightDirection, rayD);
			float proximity = lightIntensity / glm::length(mts[k].vertices[i] - light);
			float brightness = std::max(specular + incidence * proximity, brightnessThreshold);
			vBrightness.push_back(brightness);
		}

		float u = pLightRay.intersectionPoint.y;
		float v = pLightRay.intersectionPoint.z;

		pBrightness = (1.0f - u - v) * vBrightness[0] + u * vBrightness[1] + v * vBrightness[2];
	}
	else {
		pBrightness = brightnessThreshold * 0.9f;
	}

	return pBrightness;
}

float phongShading(std::vector<ModelTriangle> mts, int k, glm::vec3 rayD, float camDistance, glm::vec3 start) {
	glm::vec3 p = rayD * camDistance + start;
	glm::vec3 lightDirection = glm::normalize(p - light);
	RayTriangleIntersection lightRay = getClosestIntersection(light, lightDirection, mts);
	float pBrightness;

	if (k == lightRay.triangleIndex || mts[k].vn[0] != mts[k].vn[1]) {
		float u = lightRay.intersectionPoint.y;
		float v = lightRay.intersectionPoint.z;
		glm::vec3 normal = (1.0f - u - v) * mts[k].vn[0] + u * mts[k].vn[1] + v * mts[k].vn[2];

		float incidence = incidenceLighting(glm::normalize(normal), lightDirection);
		float specular = specularLighting(glm::normalize(normal), lightDirection, rayD);
		float proximity = lightIntensity / lightRay.distanceFromCamera;

		pBrightness = std::max(specular + incidence * proximity, brightnessThreshold);
	}
	else {
		pBrightness = brightnessThreshold * 0.9;
	}

	return pBrightness;
}

glm::vec3 refraction(glm::vec3 I, glm::vec3 N, float refIndex) {
	glm::vec3 refNormal = N;
	float etaB = 1, etaA = refIndex; // eta before and after
	float dot = glm::dot(I, refNormal);
	if (dot < 0) dot = -dot;
	else { // from medium to air
		refNormal = -N;
		std::swap(etaA, etaB);
	}

	float eta = etaB / etaA;
	float cosB = dot / (glm::length(I) * glm::length(refNormal));
	float k = 1 - eta * eta * (1 - cosB * cosB);

	if (k < 0) return glm::vec3(INT_MAX, INT_MAX, INT_MAX); // critical angle
	else return eta * I + (eta * cosB - sqrtf(k)) * refNormal;
}

Colour rayTexturing(TextureMap source, ModelTriangle t, RayTriangleIntersection lightRay) {
	std::vector<glm::vec2> texturePoints;

	for (int i = 0; i < 3; i++) {
		texturePoints.push_back(glm::vec2(t.texturePoints[i].x, t.texturePoints[i].y));
	}

	float u = lightRay.intersectionPoint.y;
	float v = lightRay.intersectionPoint.z;

	glm::vec2 intersect = (1.0f - u - v) * texturePoints[0] + u * texturePoints[1] + v * texturePoints[2];
	uint32_t c = source.pixels[int(intersect.y) * source.width + int(intersect.x)];

	uint8_t red = (c & 0x00FF0000) >> 16;
	uint8_t green = (c & 0x0000FF00) >> 8;
	uint8_t blue = (c & 0x000000FF) >> 0;
	Colour colour = Colour(red, green, blue);

	return colour;
}

uint32_t reflectiveRayTracing(std::vector<ModelTriangle> mts, glm::vec3 rayDirection, RayTriangleIntersection rayTriangle) {
	glm::vec3 reflectionPoint = rayTriangle.distanceFromCamera * rayDirection + cam_pos;
	glm::vec3 normal = rayTriangle.intersectedTriangle.normal * cam_orientation;
	glm::vec3 reflection = glm::normalize(rayDirection - 2.0f * normal * (glm::dot(rayDirection, normal)));
	RayTriangleIntersection reflectedEffect = getClosestIntersection(reflectionPoint, reflection, mts);

	if (reflectedEffect.intersectionPoint.x != INT_MAX) { // mirror shows some model
		Colour c_ref = reflectedEffect.intersectedTriangle.colour;
		Colour c = rayTriangle.intersectedTriangle.colour;
		float brightness_ref = phongShading(mts, reflectedEffect.triangleIndex, reflection, reflectedEffect.distanceFromCamera, reflectionPoint);
		float brightness = phongShading(mts, rayTriangle.triangleIndex, rayDirection, rayTriangle.distanceFromCamera, cam_pos);

		float reflectedRed = std::min(255.0f, (std::min(255.0f, c.red * brightness) * 0.3f + c_ref.red * 0.7f) * brightness_ref);
		float reflectedGreen = std::min(255.0f, (std::min(255.0f, c.green * brightness) * 0.3f + c_ref.green * 0.7f) * brightness_ref);
		float reflectedBlue = std::min(255.0f, (std::min(255.0f, c.blue * brightness) * 0.3f + c_ref.blue * 0.7f) * brightness_ref);

		uint32_t colour = (255 << 24) + (int(reflectedRed) << 16) + (int(reflectedGreen) << 8) + int(reflectedBlue);

		return colour;
	}
	else {
		Colour c = rayTriangle.intersectedTriangle.colour;
		float brightness = phongShading(mts, rayTriangle.triangleIndex, rayDirection, rayTriangle.distanceFromCamera, cam_pos);

		uint32_t colour = (255 << 24) + (int(std::min(255.0f, c.red * brightness) * 0.3f) << 16) +
			(int(std::min(255.0f, c.green * brightness) * 0.3f) << 8) + int(std::min(255.0f, c.blue * brightness) * 0.3f);

		return colour;
	}
}

uint32_t refractiveRayTracing(std::vector<ModelTriangle> mts, glm::vec3 rayDirection, RayTriangleIntersection rayTriangle) {
	glm::vec3 refractionPoint = rayTriangle.distanceFromCamera * rayDirection + cam_pos;
	glm::vec3 normal = rayTriangle.intersectedTriangle.normal * cam_orientation;
	float refractionIndex = 1.5f;
	glm::vec3 refractIn = refraction(rayDirection, normal, refractionIndex);

	std::vector<ModelTriangle> glass;
	for (int i = 12; i < 22; i++) glass.push_back(mts[i]);
	RayTriangleIntersection refractionEffect = getClosestIntersection(refractionPoint, refractIn, glass);

	if (refractionEffect.intersectionPoint.x == INT_MAX) { // refract once
		RayTriangleIntersection onceRefract = getClosestIntersection(refractionPoint, refractIn, mts);
		Colour c = onceRefract.intersectedTriangle.colour;
		float brightness = phongShading(mts, onceRefract.triangleIndex, refractIn, onceRefract.distanceFromCamera, refractionPoint);
		uint32_t colour = (255 << 24) + (int(std::min(255.0f, c.red * brightness)) << 16) +
			(int(std::min(255.0f, c.green * brightness)) << 8) + std::min(255.0f, c.blue * brightness);

		return colour;
	}
	else {
		glm::vec3 secondRefPoint = refractionEffect.distanceFromCamera * refractIn + refractionPoint;
		glm::vec3 secondNormal = refractionEffect.intersectedTriangle.normal;
		glm::vec3 refractOut = refraction(refractIn, secondNormal, refractionIndex);
		if (refractOut.x == INT_MAX) { // refract three times
			glm::vec3 secondRefIn = glm::normalize(refractIn - 2.0f * secondNormal * (glm::dot(refractIn, secondNormal)));
			RayTriangleIntersection secondRefraction = getClosestIntersection(secondRefPoint, secondRefIn, glass);
			if (secondRefraction.intersectionPoint.x != INT_MAX) {
				glm::vec3 thirdRefPoint = secondRefraction.distanceFromCamera * secondRefIn + secondRefPoint;
				glm::vec3 thirdNormal = secondRefraction.intersectedTriangle.normal;
				refractOut = refraction(secondRefIn, thirdNormal, refractionIndex);
				secondRefPoint = thirdRefPoint;
			}
			else {
				refractOut = secondRefIn;
			}
		}

		RayTriangleIntersection secondRefRay = getClosestIntersection(secondRefPoint, refractOut, mts);
		if (secondRefRay.intersectionPoint.x != INT_MAX) {
			Colour c = secondRefRay.intersectedTriangle.colour;
			float brightness = phongShading(mts, secondRefRay.triangleIndex, refractOut, secondRefRay.distanceFromCamera, secondRefPoint);
			uint32_t colour = (255 << 24) + (int(std::min(255.0f, c.red * brightness)) << 16) +
				(int(std::min(255.0f, c.green * brightness)) << 8) + std::min(255.0f, c.blue * brightness);

			return colour;
		}
	}
}

uint32_t textureRayTracing(std::vector<ModelTriangle> mts, glm::vec3 rayDirection, RayTriangleIntersection rayTriangle, std::map<std::string, Material> palette) {
	glm::vec3 p = rayDirection * rayTriangle.distanceFromCamera + cam_pos;
	glm::vec3 lightDirection = glm::normalize(p - light);
	RayTriangleIntersection lightRay = getClosestIntersection(light, lightDirection, mts);
	std::vector<glm::vec3> sourceCoords;

	float brightness = phongShading(mts, rayTriangle.triangleIndex, rayDirection, rayTriangle.distanceFromCamera, cam_pos);
	Colour c = rayTexturing(palette[mts[rayTriangle.triangleIndex].colour.name].source, rayTriangle.intersectedTriangle, lightRay);
	//brightness = 1.0f;
	uint32_t colour = (255 << 24) + (int(c.red * brightness) << 16) + (int(c.green * brightness) << 8) + c.blue * brightness;
	return colour;
}

void rayTracing(DrawingWindow& window, std::vector<ModelTriangle> mts, std::map<std::string, Material> palette) {
	float multiplier = 116.0f * 2;

	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			glm::vec3 rayVector = glm::vec3(WIDTH / 2 - x, y - HEIGHT / 2, focal * multiplier);
			glm::vec3 rayDirection = glm::normalize(-rayVector); // camera to point
			RayTriangleIntersection rayTriangle = getClosestIntersection(cam_pos, rayDirection, mts);

			if (rayTriangle.intersectionPoint.x != INT_MAX) { // ray hits the model
				if (rayTriangle.triangleIndex == 26 || rayTriangle.triangleIndex == 31) { // reflective surface
					uint32_t colour = reflectiveRayTracing(mts, rayDirection, rayTriangle);
					window.setPixelColour(x, y, colour);
					window.renderFrame();
				}
				else if (rayTriangle.triangleIndex >= 12 && rayTriangle.triangleIndex <= 21) { // refractive material
					uint32_t colour = refractiveRayTracing(mts, rayDirection, rayTriangle);
					window.setPixelColour(x, y, colour);
					window.renderFrame();
				}
				else if (palette[mts[rayTriangle.triangleIndex].colour.name].texturePresent == true) { // raytracing texturing
					uint32_t colour = textureRayTracing(mts, rayDirection, rayTriangle, palette);
					window.setPixelColour(x, y, colour);
					window.renderFrame();
				}
				else {
					Colour c = rayTriangle.intersectedTriangle.colour;
					float brightness = phongShading(mts, rayTriangle.triangleIndex, rayDirection, rayTriangle.distanceFromCamera, cam_pos);
					uint32_t colour = (255 << 24) + (int(std::min(255.0f, c.red * brightness)) << 16) +
						(int(std::min(255.0f, c.green * brightness)) << 8) + int(std::min(255.0f, c.blue * brightness));
					window.setPixelColour(x, y, colour);
					window.renderFrame();
				}
			}
		}
	}

	std::cout << "rayTracing done" << std::endl;
}

void rasterising(DrawingWindow& window, std::vector<ModelTriangle> mts, std::map<std::string, Material> palette) {
	float multiplier = 116.0f * 2;
	for (int i = 0; i < mts.size(); i++) {
		CanvasPoint v0;
		CanvasPoint v1;
		CanvasPoint v2;
		std::vector<CanvasPoint> t_plane = { v0, v1, v2 };
		for (int k = 0; k < 3; k++) {
			CanvasPoint v_model = CanvasPoint();
			glm::vec3 coor_m = glm::vec3(mts[i].vertices[k].x, mts[i].vertices[k].y, mts[i].vertices[k].z);
			//coor_m = (cam_pos - coor_m) * cam_orientation;
			coor_m = cam_pos - coor_m * cam_orientation; // rotate the source
			float z_depth = 1 / coor_m.z;
			v_model.x = -(multiplier * focal * coor_m.x / coor_m.z) + WIDTH / 2;
			v_model.y = multiplier * focal * coor_m.y / coor_m.z + HEIGHT / 2;
			t_plane[k] = CanvasPoint(v_model.x, v_model.y, z_depth);
		}
		CanvasTriangle t = CanvasTriangle(t_plane[0], t_plane[1], t_plane[2]);
		t.vertices[0].texturePoint = mts[i].texturePoints[0];
		t.vertices[1].texturePoint = mts[i].texturePoints[1];
		t.vertices[2].texturePoint = mts[i].texturePoints[2];

		if (palette[mts[i].colour.name].texturePresent == false) drawFilledTriangle(window, t, mts[i]);
		else drawTextureTriangle(window, palette[mts[i].colour.name].source, t);
	}
}

int main(int argc, char* argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	float scale = 1.2f;
	std::map<std::string, Material> palette = loadMtl("textured-cornell-box.mtl");
	std::vector<ModelTriangle> mts1 = loadOBJ("cornell-box.obj", palette, scale);
	std::vector<ModelTriangle> mts2 = loadOBJ("sphere.obj", palette, scale);
	//std::vector<ModelTriangle> mts3 = loadOBJ("complexSphere.obj", palette, scale);
	for (int i = 0; i < mts2.size(); i++) mts1.push_back(mts2[i]);
	std::vector<ModelTriangle> mts;
	//for (int i = 0; i < 20; i++) mts.push_back(mts1[i]);
	//for (int i = 8; i < mts1.size(); i++) mts.push_back(mts1[i]);
	//uint32_t test = (255 << 24) + (255 << 16) + (255 << 8) + 255;
	//window.setPixelColour(-10, 20, test);

	/*uint32_t c = (78 << 24) + (193 << 16) + (255 << 8) + 136;
	std::cout << "c: " << c << std::endl;

	glm::vec3 vertex = glm::vec3(13.0f, -39.0f, 21.0f);
	vertex = cam_pos - vertex;
	float x = focal * vertex.x / vertex.z + WIDTH / 2;
	float y = focal * vertex.y / vertex.z + HEIGHT / 2;
	std::cout << "x: " << x << std::endl;
	std::cout << "y: " << y << std::endl;

	glm::vec3 v0 = glm::vec3(304.0f, 61.0f, 21.0f);
	glm::vec3 v1 = glm::vec3(498.0f, 608.0f, 21.0f);
	glm::vec3 v2 = glm::vec3(57.0f, 404.0f, 21.0f);

	glm::vec3 p = glm::vec3(286.0f, 357.0f, 21.0f);

	glm::vec3 vt0 = glm::vec3(317.0f, 33.0f, 21.0f);
	glm::vec3 vt1 = glm::vec3(515.0f, 588.0f, 21.0f);
	glm::vec3 vt2 = glm::vec3(50.0f, 402.0f, 21.0f);

	float ap = glm::length(glm::cross((v0 - v1), (v0 - v2)));
	float at = glm::length(glm::cross((vt0 - vt1), (vt0 - vt2)));
	float u = glm::length(glm::cross((p - v0), (p - v2)))/ap;
	float v = glm::length(glm::cross((p - v1), (p - v2)))/ap;
	float w = glm::length(glm::cross((p - v0), (p - v1)))/ap;

	glm::vec3 vt = u*vt1 + v*vt0 + w*vt2;
	std::cout << "x: " << vt.x << std::endl;
	std::cout << "y: " << vt.y << std::endl;*/

	/*glm::vec3 cam_position = glm::vec3(61.0f, 18.0f, -17.0f);
	glm::vec3 object = glm::vec3(-61.0f, 54.0f, -22.0f);
	glm::vec3 z = glm::normalize(cam_pos - object);
	glm::vec3 x = glm::cross({ 0.0f, 1.0f, 0.0f }, z);
	glm::vec3 y = glm::cross(z, x);

	glm::mat3 orientMatrix = { x, y, z };*/

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		/*window.clearPixels();
		for (size_t y = 0; y < window.height; y++) { for (size_t x = 0; x < window.width; x++) depthMap[y][x] = 0.0f; }
		rasterising(window, mts1, palette);*/
		rayTracing(window, mts1, palette);
		//rotateSource(rotateY(0.2));
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		//window.renderFrame();
	}
}
