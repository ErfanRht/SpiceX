#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>
#ifdef _WIN32
#include <SDL2/SDL2_gfx.h>
#else
#include <SDL2/SDL2_gfxPrimitives.h>
#endif
#include <bits/stdc++.h>
#include <cmath>

using namespace std;

const int SCREEN_WIDTH = 1280;
const int SCREEN_HEIGHT = 720;
const int TOP_BAR_HEIGHT = 40;
const int GRID_SPACING = 20;
const int COMPONENT_DEFAULT_LENGTH = GRID_SPACING * 4;

const string BASE_PATH = "C:/Users/Erfan/Dev/Cpp/sutSpice_phase2/";
const string ASSET_PATH = BASE_PATH + "assets/";
const string SCHEMATICS_PATH = BASE_PATH + "schematics/";

bool isFileMenuOpen = false;
bool isComponentsMenuOpen = false;
bool isSimulateMenuOpen = false;

enum class InteractionMode {
    NONE,
    WIRING,
    PLACE_COMPONENT,
    DRAGGING_COMPONENT,
    DELETE_ITEM,
    EDITING_COMPONENT_VALUE,
    PLACE_LABEL,
    EDITING_LABEL_TEXT,
    PLACE_GND_LABEL,
    DIALOG_ACTIVE
};
InteractionMode currentInteractionMode = InteractionMode::NONE;
bool isDrawingWire = false;
SDL_Point firstInteractionPoint = {-1, -1};

enum class ComponentType {
    NONE, RESISTOR, CAPACITOR, INDUCTOR, DIODE, DC_VOLTAGE_SOURCE, DC_CURRENT_SOURCE, AC_VOLTAGE_SOURCE
};
ComponentType selectedComponentType = ComponentType::NONE;
int placementRotation = 0;

struct Wire {
    SDL_Point start;
    SDL_Point end;
};

struct Component {
    ComponentType type;
    SDL_Point node1;
    SDL_Point node2;
    string id;
    int rotation_angle = 0;
    bool is_selected = false;
    string value = "R";
    SDL_Rect value_rect;
    SDL_Rect id_rect;
};

struct NodeLabel {
    string text;
    SDL_Point position;
    SDL_Rect text_rect;
    int rotation_angle = 0;
};

vector<Wire> wires;
vector<Component> components;
vector<NodeLabel> labels;
Component* selectedComponent = nullptr;
Component* editingComponent = nullptr;
NodeLabel* editingLabel = nullptr;
string textInputBuffer = "";
int drag_offset_x = 0;
int drag_offset_y = 0;

string currentSchematicFileName = "";
bool schematicModified = false;

map<ComponentType, int> componentNameCounters;

enum class DialogType { NONE, SAVE_AS, OPEN, TRANSIENT_ANALYSIS, DC_SWEEP_ANALYSIS };
DialogType activeDialogType = DialogType::NONE;

struct DialogField {
    string label;
    string buffer;
    SDL_Rect input_rect;
};
vector<DialogField> dialogFields;
string dialogTitle;
int activeDialogFieldIndex = -1;

class Button {
public:
    Button(string text, int x, int y, int w, int h, function<void()> onClick) {
        m_text = text;
        m_position = {x, y, w, h};
        m_onClick = onClick;
        m_is_hovered = false;
    }

    SDL_Rect m_position;
    function<void()> m_onClick;
    bool m_is_hovered;

    void handle_event(SDL_Event* e) {
        if (e->type == SDL_MOUSEBUTTONDOWN) {
            int x, y;
            SDL_GetMouseState(&x, &y);
            if (x > m_position.x && x < m_position.x + m_position.w &&
                y > m_position.y && y < m_position.y + m_position.h) {
                m_onClick();
            }
        } else if (e->type == SDL_MOUSEMOTION) {
            int x, y;
            SDL_GetMouseState(&x, &y);
            m_is_hovered = (x > m_position.x && x < m_position.x + m_position.w &&
                            y > m_position.y && y < m_position.y + m_position.h);
        }
    }

    void render(SDL_Renderer* renderer, TTF_Font* font) {
        if (m_is_hovered) {
            SDL_SetRenderDrawColor(renderer, 0x66, 0x66, 0x66, 0xFF);
        } else {
            SDL_SetRenderDrawColor(renderer, 0x44, 0x44, 0x44, 0xFF);
        }
        SDL_RenderFillRect(renderer, &m_position);

        if (font && !m_text.empty()) {
            SDL_Color textColor = { 0xFF, 0xFF, 0xFF, 0xFF };
            SDL_Surface* textSurface = TTF_RenderText_Solid(font, m_text.c_str(), textColor);
            if (textSurface) {
                SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                SDL_Rect renderQuad = { m_position.x + (m_position.w - textSurface->w) / 2, m_position.y + (m_position.h - textSurface->h) / 2, textSurface->w, textSurface->h };
                SDL_RenderCopy(renderer, textTexture, nullptr, &renderQuad);
                SDL_DestroyTexture(textTexture);
                SDL_FreeSurface(textSurface);
            }
        }
    }
private:
    string m_text;
};

struct ComponentMenuItem {
    string name;
    ComponentType type;
    SDL_Texture* iconTexture = nullptr;
    Button button;
};

SDL_Point snap_to_grid(int x, int y) {
    int snappedX = round((float)x / GRID_SPACING) * GRID_SPACING;
    int snappedY = round((float)y / GRID_SPACING) * GRID_SPACING;
    if (snappedY < TOP_BAR_HEIGHT) {
        snappedY = TOP_BAR_HEIGHT;
    }
    return { snappedX, snappedY };
}

void draw_grid(SDL_Renderer* renderer) {
    SDL_SetRenderDrawColor(renderer, 0x40, 0x40, 0x40, 0xFF);
    for (int x = 0; x < SCREEN_WIDTH; x += GRID_SPACING) {
        for (int y = TOP_BAR_HEIGHT; y < SCREEN_HEIGHT; y += GRID_SPACING) {
            SDL_RenderDrawPoint(renderer, x, y);
        }
    }
}

SDL_Texture* loadAndProcessTexture(const string& path, SDL_Renderer* renderer) {
    SDL_Surface* originalSurface = IMG_Load(path.c_str());
    if (!originalSurface) { cerr << "Failed to load image " << path << "! SDL_image Error: " << IMG_GetError() << endl; return nullptr; }
    SDL_LockSurface(originalSurface);
    auto pixels = static_cast<uint32_t*>(originalSurface->pixels);
    int pixelCount = originalSurface->w * originalSurface->h;
    SDL_PixelFormat* format = originalSurface->format;
    uint32_t white = SDL_MapRGB(format, 255, 255, 255);
    for (int i = 0; i < pixelCount; ++i) {
        uint8_t r, g, b, a;
        SDL_GetRGBA(pixels[i], format, &r, &g, &b, &a);
        if (a > 50) pixels[i] = white;
    }
    SDL_UnlockSurface(originalSurface);
    SDL_Texture* finalTexture = SDL_CreateTextureFromSurface(renderer, originalSurface);
    SDL_FreeSurface(originalSurface);
    SDL_SetTextureBlendMode(finalTexture, SDL_BLENDMODE_BLEND);
    return finalTexture;
}

void get_component_defaults(ComponentType type, string& value) {
    switch (type) {
        case ComponentType::RESISTOR:            value = "1k"; break;
        case ComponentType::CAPACITOR:           value = "1u"; break;
        case ComponentType::INDUCTOR:            value = "1m"; break;
        case ComponentType::DIODE:               value = "1N4148"; break;
        case ComponentType::DC_VOLTAGE_SOURCE:   value = "5"; break;
        case ComponentType::DC_CURRENT_SOURCE:   value = "1m"; break;
        case ComponentType::AC_VOLTAGE_SOURCE:   value = "1"; break;
        default:                                 value = "?"; break;
    }
}

string generate_component_name(ComponentType type) {
    char prefix = '?';
    switch (type) {
        case ComponentType::RESISTOR:            prefix = 'R'; break;
        case ComponentType::CAPACITOR:           prefix = 'C'; break;
        case ComponentType::INDUCTOR:            prefix = 'L'; break;
        case ComponentType::DIODE:               prefix = 'D'; break;
        case ComponentType::DC_VOLTAGE_SOURCE:   prefix = 'V'; break;
        case ComponentType::DC_CURRENT_SOURCE:   prefix = 'I'; break;
        case ComponentType::AC_VOLTAGE_SOURCE:   prefix = 'V'; break;
        default: break;
    }
    componentNameCounters[type]++;
    return prefix + to_string(componentNameCounters[type]);
}


string component_type_to_string(ComponentType type) {
    switch (type) {
        case ComponentType::RESISTOR: return "RESISTOR";
        case ComponentType::CAPACITOR: return "CAPACITOR";
        case ComponentType::INDUCTOR: return "INDUCTOR";
        case ComponentType::DIODE: return "DIODE";
        case ComponentType::DC_VOLTAGE_SOURCE: return "DC_VOLTAGE_SOURCE";
        case ComponentType::DC_CURRENT_SOURCE: return "DC_CURRENT_SOURCE";
        case ComponentType::AC_VOLTAGE_SOURCE: return "AC_VOLTAGE_SOURCE";
        default: return "UNKNOWN";
    }
}

ComponentType string_to_component_type(const string& str) {
    if (str == "RESISTOR") return ComponentType::RESISTOR;
    if (str == "CAPACITOR") return ComponentType::CAPACITOR;
    if (str == "INDUCTOR") return ComponentType::INDUCTOR;
    if (str == "DIODE") return ComponentType::DIODE;
    if (str == "DC_VOLTAGE_SOURCE") return ComponentType::DC_VOLTAGE_SOURCE;
    if (str == "DC_CURRENT_SOURCE") return ComponentType::DC_CURRENT_SOURCE;
    if (str == "AC_VOLTAGE_SOURCE") return ComponentType::AC_VOLTAGE_SOURCE;
    return ComponentType::NONE;
}

void update_name_counters_from_id(ComponentType type, const string& id) {
    if (id.empty() || !isdigit(id.back())) return;
    size_t first_digit = id.find_first_of("0123456789");
    if (first_digit == string::npos) return;
    int num = stoi(id.substr(first_digit));
    if (num > componentNameCounters[type]) {
        componentNameCounters[type] = num;
    }
}

void save_schematic(const string& filename) {
    ofstream outFile(filename);
    if (!outFile.is_open()) { cerr << "Error: Could not open file for saving: " << filename << endl; return; }
    for (const auto& comp : components) {
        outFile << "COMPONENT," << component_type_to_string(comp.type) << "," << comp.id << ","
                << comp.node1.x << "," << comp.node1.y << "," << comp.node2.x << "," << comp.node2.y << ","
                << comp.rotation_angle << "," << comp.value << endl;
    }
    for (const auto& wire : wires) {
        outFile << "WIRE," << wire.start.x << "," << wire.start.y << "," << wire.end.x << "," << wire.end.y << endl;
    }
    for (const auto& label : labels) {
        outFile << "LABEL," << label.text << "," << label.position.x << "," << label.position.y << "," << label.rotation_angle << endl;
    }
    outFile.close();
    currentSchematicFileName = filename;
    schematicModified = false;
    cout << "Schematic saved to: " << filename << endl;
}

void load_schematic(const string& filename) {
    ifstream inFile(filename);
    if (!inFile.is_open()) { cerr << "Error: Could not open file for loading: " << filename << endl; return; }
    components.clear();
    wires.clear();
    labels.clear();
    componentNameCounters.clear();
    string line;
    while (getline(inFile, line)) {
        stringstream ss(line);
        string type_str;
        getline(ss, type_str, ',');
        if (type_str == "COMPONENT") {
            string comp_type_str, id_str, x1_str, y1_str, x2_str, y2_str, rot_str, val_str;
            getline(ss, comp_type_str, ','); getline(ss, id_str, ','); getline(ss, x1_str, ','); getline(ss, y1_str, ',');
            getline(ss, x2_str, ','); getline(ss, y2_str, ','); getline(ss, rot_str, ','); getline(ss, val_str);
            ComponentType type = string_to_component_type(comp_type_str);
            components.push_back({type, {stoi(x1_str), stoi(y1_str)}, {stoi(x2_str), stoi(y2_str)}, id_str, stoi(rot_str), false, val_str});
            update_name_counters_from_id(type, id_str);
        } else if (type_str == "WIRE") {
            string x1_str, y1_str, x2_str, y2_str;
            getline(ss, x1_str, ','); getline(ss, y1_str, ','); getline(ss, x2_str, ','); getline(ss, y2_str);
            wires.push_back({{stoi(x1_str), stoi(y1_str)}, {stoi(x2_str), stoi(y2_str)}});
        } else if (type_str == "LABEL") {
            string text, x_str, y_str, rot_str;
            getline(ss, text, ','); getline(ss, x_str, ','); getline(ss, y_str, ','); getline(ss, rot_str);
            labels.push_back({text, {stoi(x_str), stoi(y_str)}, {}, stoi(rot_str)});
        }
    }
    inFile.close();
    currentSchematicFileName = filename;
    schematicModified = false;
    cout << "Schematic loaded from: " << filename << endl;
}

void draw_schematic_elements(SDL_Renderer* renderer, TTF_Font* valueFont, const vector<ComponentMenuItem>& componentMenuIcons) {
    for (const auto& wire : wires) {
        thickLineRGBA(renderer, wire.start.x, wire.start.y, wire.end.x, wire.end.y, 3, 0xFF, 0xFF, 0xFF, 0xFF);
        filledCircleRGBA(renderer, wire.start.x, wire.start.y, 0, 220, 20, 60, 0xFF);
        filledCircleRGBA(renderer, wire.end.x, wire.end.y, 0, 220, 20, 60, 0xFF);
    }
    for (auto& comp : components) {
        SDL_Texture* iconTexture = nullptr;
        for (const auto& item : componentMenuIcons) if (item.type == comp.type) { iconTexture = item.iconTexture; break; }

        int center_x = (comp.node1.x + comp.node2.x) / 2;
        int center_y = (comp.node1.y + comp.node2.y) / 2;

        if (iconTexture) {
            int icon_size = COMPONENT_DEFAULT_LENGTH;
            int icon_half_size = icon_size / 2;
            SDL_Rect destRect = { center_x - icon_half_size, center_y - icon_half_size, icon_size, icon_size };
            SDL_Point rotationCenter = { icon_half_size, icon_half_size };
            SDL_RenderCopyEx(renderer, iconTexture, nullptr, &destRect, comp.rotation_angle, &rotationCenter, SDL_FLIP_NONE);
        }
        filledCircleRGBA(renderer, comp.node1.x, comp.node1.y, 7, 220, 20, 60, 0xFF);
        filledCircleRGBA(renderer, comp.node2.x, comp.node2.y, 7, 220, 20, 60, 0xFF);

        SDL_Color idColor = {0xAA, 0xAA, 0xFF, 0xFF};
        SDL_Surface* idSurface = TTF_RenderText_Solid(valueFont, comp.id.c_str(), idColor);
        if (idSurface) {
            SDL_Texture* idTexture = SDL_CreateTextureFromSurface(renderer, idSurface);
            comp.id_rect = {center_x - idSurface->w / 2, center_y - 50, idSurface->w, idSurface->h};
            SDL_RenderCopy(renderer, idTexture, nullptr, &comp.id_rect);
            SDL_FreeSurface(idSurface);
            SDL_DestroyTexture(idTexture);
        }

        string text_to_render = (editingComponent == &comp) ? textInputBuffer : comp.value;
        SDL_Color textColor = { 0xFF, 0xFF, 0xFF, 0xFF };
        SDL_Surface* textSurface = TTF_RenderText_Solid(valueFont, text_to_render.c_str(), textColor);
        if (textSurface) {
            SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
            comp.value_rect = {center_x - textSurface->w / 2, center_y + 35, textSurface->w, textSurface->h};
            if (editingComponent == &comp) {
                SDL_SetRenderDrawColor(renderer, 0x00, 0x00, 0x55, 0xFF);
                SDL_RenderFillRect(renderer, &comp.value_rect);
            }
            SDL_RenderCopy(renderer, textTexture, nullptr, &comp.value_rect);
            SDL_FreeSurface(textSurface);
            SDL_DestroyTexture(textTexture);
        }
    }
    for(auto& label : labels) {
        string text_to_render = (editingLabel == &label) ? textInputBuffer : label.text;
        SDL_Color textColor = (label.text == "GND") ? SDL_Color{100, 255, 100, 255} : SDL_Color{0x99, 0xFF, 0xFF, 0xFF};
        SDL_Surface* textSurface = TTF_RenderText_Solid(valueFont, text_to_render.c_str(), textColor);
        if(textSurface) {
            SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
            int w = textSurface->w;
            int h = textSurface->h;
            int offsetX = 10;
            int offsetY = 10;

            if (label.rotation_angle == 0) { label.text_rect = {label.position.x + offsetX, label.position.y - h - (offsetY/2), w, h}; }
            else if (label.rotation_angle == 90) { label.text_rect = {label.position.x - w - offsetX, label.position.y - h/2, w, h}; }
            else if (label.rotation_angle == 180) { label.text_rect = {label.position.x - w/2, label.position.y + offsetY, w, h}; }
            else { label.text_rect = {label.position.x + offsetX, label.position.y - h/2, w, h}; }

            if (editingLabel == &label) {
                SDL_SetRenderDrawColor(renderer, 0x00, 0x33, 0x33, 0xFF);
                SDL_RenderFillRect(renderer, &label.text_rect);
            }
            SDL_RenderCopy(renderer, textTexture, nullptr, &label.text_rect);
            SDL_FreeSurface(textSurface);
            SDL_DestroyTexture(textTexture);
        }
    }
}

double dist_to_segment_sq(SDL_Point p, SDL_Point v, SDL_Point w) {
    double l2 = pow(v.x - w.x, 2) + pow(v.y - w.y, 2);
    if (l2 == 0.0) return pow(p.x - v.x, 2) + pow(p.y - v.y, 2);
    double t = max(0.0, min(1.0, ((p.x - v.x) * (w.x - v.x) + (p.y - v.y) * (w.y - v.y)) / l2));
    return pow(p.x - (v.x + t * (w.x - v.x)), 2) + pow(p.y - (v.y + t * (w.y - v.y)), 2);
}

void close_dialog() {
    activeDialogType = DialogType::NONE;
    dialogFields.clear();
    dialogTitle = "";
    activeDialogFieldIndex = -1;
    currentInteractionMode = InteractionMode::NONE;
    SDL_StopTextInput();
}

void setup_dialog(DialogType type, const string& title, const vector<string>& labels) {
    activeDialogType = type;
    dialogTitle = title;
    dialogFields.clear();
    for (const auto& label : labels) {
        dialogFields.push_back({label, ""});
    }
    activeDialogFieldIndex = 0;
    currentInteractionMode = InteractionMode::DIALOG_ACTIVE;
    SDL_StartTextInput();
}

void handle_dialog_ok() {
    if (activeDialogType == DialogType::SAVE_AS) {
        if (!dialogFields.empty() && !dialogFields[0].buffer.empty()) {
            save_schematic(SCHEMATICS_PATH + dialogFields[0].buffer);
        }
    } else if (activeDialogType == DialogType::OPEN) {
        if (!dialogFields.empty() && !dialogFields[0].buffer.empty()) {
            load_schematic(SCHEMATICS_PATH + dialogFields[0].buffer);
        }
    } else if (activeDialogType == DialogType::TRANSIENT_ANALYSIS) {
        cout << "--- TRANSIENT ANALYSIS PARAMETERS ---" << endl;
        cout << "Tstep: " << dialogFields[0].buffer << endl;
        cout << "Tstop: " << dialogFields[1].buffer << endl;
        cout << "Wanted Value: " << dialogFields[2].buffer << endl;
        cout << "------------------------------------" << endl;
    } else if (activeDialogType == DialogType::DC_SWEEP_ANALYSIS) {
        cout << "--- DC SWEEP ANALYSIS PARAMETERS ---" << endl;
        cout << "Source Name: " << dialogFields[0].buffer << endl;
        cout << "Start Value: " << dialogFields[1].buffer << endl;
        cout << "End Value: " << dialogFields[2].buffer << endl;
        cout << "Increment: " << dialogFields[3].buffer << endl;
        cout << "------------------------------------" << endl;
    }
    close_dialog();
}

void render_dialog(SDL_Renderer* renderer, TTF_Font* font) {
    const int dialog_w = 450;
    const int dialog_h = 100 + dialogFields.size() * 40 + 50;
    SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };

    SDL_SetRenderDrawColor(renderer, 0x55, 0x55, 0x55, 0xFF);
    SDL_RenderFillRect(renderer, &panelRect);
    SDL_SetRenderDrawColor(renderer, 0x88, 0x88, 0x88, 0xFF);
    SDL_RenderDrawRect(renderer, &panelRect);

    SDL_Color textColor = { 0xFF, 0xFF, 0xFF, 0xFF };
    SDL_Surface* titleSurface = TTF_RenderText_Solid(font, dialogTitle.c_str(), textColor);
    if (titleSurface) {
        SDL_Texture* titleTexture = SDL_CreateTextureFromSurface(renderer, titleSurface);
        SDL_Rect titleRect = { panelRect.x + (panelRect.w - titleSurface->w) / 2, panelRect.y + 20, titleSurface->w, titleSurface->h };
        SDL_RenderCopy(renderer, titleTexture, nullptr, &titleRect);
        SDL_FreeSurface(titleSurface);
        SDL_DestroyTexture(titleTexture);
    }

    int current_y = panelRect.y + 60;
    for (size_t i = 0; i < dialogFields.size(); ++i) {
        SDL_Surface* labelSurface = TTF_RenderText_Solid(font, dialogFields[i].label.c_str(), textColor);
        if (labelSurface) {
            SDL_Texture* labelTexture = SDL_CreateTextureFromSurface(renderer, labelSurface);
            SDL_Rect labelRect = { panelRect.x + 20, current_y + (30 - labelSurface->h)/2, labelSurface->w, labelSurface->h };
            SDL_RenderCopy(renderer, labelTexture, nullptr, &labelRect);
            SDL_FreeSurface(labelSurface);
            SDL_DestroyTexture(labelTexture);
        }

        SDL_Rect inputRect = { panelRect.x + 150, current_y, panelRect.w - 170, 30 };
        dialogFields[i].input_rect = inputRect;
        if (activeDialogFieldIndex == i) {
            SDL_SetRenderDrawColor(renderer, 0x00, 0x00, 0x44, 0xFF);
        } else {
            SDL_SetRenderDrawColor(renderer, 0x33, 0x33, 0x33, 0xFF);
        }
        SDL_RenderFillRect(renderer, &inputRect);
        SDL_SetRenderDrawColor(renderer, 0xAAAAAA, 0xAAAAAA, 0xAAAAAA, 0xFF);
        SDL_RenderDrawRect(renderer, &inputRect);

        if (!dialogFields[i].buffer.empty()) {
            SDL_Surface* textSurface = TTF_RenderText_Solid(font, dialogFields[i].buffer.c_str(), textColor);
            if (textSurface) {
                SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                SDL_Rect textRect = { inputRect.x + 5, inputRect.y + (inputRect.h - textSurface->h) / 2, textSurface->w, textSurface->h };
                SDL_RenderCopy(renderer, textTexture, nullptr, &textRect);
                SDL_FreeSurface(textSurface);
                SDL_DestroyTexture(textTexture);
            }
        }
        current_y += 40;
    }

    Button okButton("OK", panelRect.x + panelRect.w - 190, panelRect.y + panelRect.h - 45, 80, 30, handle_dialog_ok);
    Button cancelButton("Cancel", panelRect.x + panelRect.w - 100, panelRect.y + panelRect.h - 45, 80, 30, close_dialog);

    int mx, my;
    SDL_GetMouseState(&mx, &my);
    SDL_Point mousePoint = {mx, my};
    okButton.m_is_hovered = SDL_PointInRect(&mousePoint, &okButton.m_position);
    cancelButton.m_is_hovered = SDL_PointInRect(&mousePoint, &cancelButton.m_position);

    okButton.render(renderer, font);
    cancelButton.render(renderer, font);
}


int main(int argc, char* args[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0 || TTF_Init() == -1) { cerr << "Core SDL init failed." << endl; return -1; }
    if (!(IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG)) { cerr << "SDL_image init failed." << endl; return -1; }
    SDL_Window* window = SDL_CreateWindow("SUTSpice | Phase 2 (Updated)", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    TTF_Font* uiFont = TTF_OpenFont(R"(C:\Windows\Fonts\arial.ttf)", 16);
    TTF_Font* valueFont = TTF_OpenFont(R"(C:\Windows\Fonts\arialbd.ttf)", 18);
    if (!window || !renderer || !uiFont || !valueFont) { cerr << "Window or font creation failed." << endl; return -1; }

    vector<Button> topBarButtons;
    topBarButtons.emplace_back("File", 10, 5, 80, 30, [](){ isFileMenuOpen = !isFileMenuOpen; isComponentsMenuOpen = false; isSimulateMenuOpen = false; currentInteractionMode = InteractionMode::NONE; });
    topBarButtons.emplace_back("Simulate", 100, 5, 100, 30, [](){ isSimulateMenuOpen = !isSimulateMenuOpen; isFileMenuOpen = false; isComponentsMenuOpen = false; currentInteractionMode = InteractionMode::NONE; });
    topBarButtons.emplace_back("Add Component", 210, 5, 150, 30, [](){ isComponentsMenuOpen = !isComponentsMenuOpen; isFileMenuOpen = false; isSimulateMenuOpen = false; currentInteractionMode = InteractionMode::NONE; });
    topBarButtons.emplace_back("Add Wire", 370, 5, 100, 30, [](){ currentInteractionMode = InteractionMode::WIRING; cout << "Mode: Add Wire" << endl; });
    topBarButtons.emplace_back("Add GND", 480, 5, 100, 30, [](){ currentInteractionMode = InteractionMode::PLACE_GND_LABEL; cout << "Mode: Add GND" << endl; });
    topBarButtons.emplace_back("Add Label", 590, 5, 100, 30, [](){ currentInteractionMode = InteractionMode::PLACE_LABEL; cout << "Mode: Add Label" << endl; });
    topBarButtons.emplace_back("Delete", 700, 5, 100, 30, [](){ currentInteractionMode = InteractionMode::DELETE_ITEM; cout << "Mode: Delete Item" << endl; });

    vector<Button> fileMenuButtons;
    fileMenuButtons.emplace_back("New Schematic", 10, 40, 180, 30, [](){ components.clear(); wires.clear(); labels.clear(); componentNameCounters.clear(); currentSchematicFileName = ""; schematicModified = false; isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Open...", 10, 70, 180, 30, [](){ setup_dialog(DialogType::OPEN, "Open Schematic", {"Filename"}); isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Save", 10, 100, 180, 30, [](){ if (currentSchematicFileName.empty()) { setup_dialog(DialogType::SAVE_AS, "Save Schematic As", {"Filename"}); } else { save_schematic(currentSchematicFileName); } isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Save As...", 10, 130, 180, 30, [](){ setup_dialog(DialogType::SAVE_AS, "Save Schematic As", {"Filename"}); isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Clear", 10, 160, 180, 30, [](){ components.clear(); wires.clear(); labels.clear(); componentNameCounters.clear(); schematicModified = true; isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Exit", 10, 190, 180, 30, [](){  });

    vector<Button> simulateMenuButtons;
    simulateMenuButtons.emplace_back("DC Sweep Analysis", 100, 40, 180, 30, [](){ setup_dialog(DialogType::DC_SWEEP_ANALYSIS, "DC Sweep Analysis", {"Source Name", "Start Value", "End Value", "Increment"}); isSimulateMenuOpen = false; });
    simulateMenuButtons.emplace_back("Transient Analysis", 100, 70, 180, 30, [](){ setup_dialog(DialogType::TRANSIENT_ANALYSIS, "Transient Analysis", {"Tstep", "Tstop", "Wanted Value"}); isSimulateMenuOpen = false; });

    vector<ComponentMenuItem> componentMenuItems;
    const int COMP_ITEM_W = 150, COMP_ITEM_H = 100, COMP_ITEM_PAD = 10;
    auto selectComp = [](ComponentType t){ selectedComponentType = t; currentInteractionMode = InteractionMode::PLACE_COMPONENT; isComponentsMenuOpen = false; placementRotation = 0; };
    componentMenuItems.push_back({"Resistor", ComponentType::RESISTOR, loadAndProcessTexture(ASSET_PATH + "resistor.png", renderer), Button("", 210, 40, COMP_ITEM_W, COMP_ITEM_H, [=](){ selectComp(ComponentType::RESISTOR); })});
    componentMenuItems.push_back({"Capacitor", ComponentType::CAPACITOR, loadAndProcessTexture(ASSET_PATH + "capacitor.png", renderer), Button("", 210 + (COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=](){ selectComp(ComponentType::CAPACITOR); })});
    componentMenuItems.push_back({"Inductor", ComponentType::INDUCTOR, loadAndProcessTexture(ASSET_PATH + "inductor.png", renderer), Button("", 210 + 2*(COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=](){ selectComp(ComponentType::INDUCTOR); })});
    componentMenuItems.push_back({"Diode", ComponentType::DIODE, loadAndProcessTexture(ASSET_PATH + "diode.png", renderer), Button("", 210 + 3*(COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=](){ selectComp(ComponentType::DIODE); })});
    componentMenuItems.push_back({"DC Voltage Source", ComponentType::DC_VOLTAGE_SOURCE, loadAndProcessTexture(ASSET_PATH + "voltage_source.png", renderer), Button("", 210, 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=](){ selectComp(ComponentType::DC_VOLTAGE_SOURCE); })});
    componentMenuItems.push_back({"DC Current Source", ComponentType::DC_CURRENT_SOURCE, loadAndProcessTexture(ASSET_PATH + "current_source.png", renderer), Button("", 210 + (COMP_ITEM_W + COMP_ITEM_PAD), 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=](){ selectComp(ComponentType::DC_CURRENT_SOURCE); })});
    componentMenuItems.push_back({"AC Voltage Source", ComponentType::AC_VOLTAGE_SOURCE, loadAndProcessTexture(ASSET_PATH + "ac_voltage_source.png", renderer), Button("", 210 + 2*(COMP_ITEM_W + COMP_ITEM_PAD), 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=](){ selectComp(ComponentType::AC_VOLTAGE_SOURCE); })});

    bool quit = false;
    SDL_Event e;
    while (!quit) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) { quit = true; }

            if (currentInteractionMode == InteractionMode::DIALOG_ACTIVE) {
                if (e.type == SDL_KEYDOWN) {
                    if (e.key.keysym.sym == SDLK_RETURN) { handle_dialog_ok(); }
                    else if (e.key.keysym.sym == SDLK_ESCAPE) { close_dialog(); }
                    else if (e.key.keysym.sym == SDLK_TAB) {
                        if(activeDialogFieldIndex != -1) {
                            activeDialogFieldIndex = (activeDialogFieldIndex + 1) % dialogFields.size();
                        }
                    }
                    else if (e.key.keysym.sym == SDLK_BACKSPACE && activeDialogFieldIndex != -1 && !dialogFields[activeDialogFieldIndex].buffer.empty()) {
                        dialogFields[activeDialogFieldIndex].buffer.pop_back();
                    }
                } else if (e.type == SDL_TEXTINPUT && activeDialogFieldIndex != -1) {
                    dialogFields[activeDialogFieldIndex].buffer += e.text.text;
                } else if (e.type == SDL_MOUSEBUTTONDOWN) {
                    int mx, my;
                    SDL_GetMouseState(&mx, &my);
                    SDL_Point mousePoint = {mx, my};
                    bool fieldClicked = false;
                    for(size_t i = 0; i < dialogFields.size(); ++i) {
                        if(SDL_PointInRect(&mousePoint, &dialogFields[i].input_rect)) {
                            activeDialogFieldIndex = i;
                            fieldClicked = true;
                            break;
                        }
                    }
                    const int dialog_w = 450;
                    const int dialog_h = 100 + dialogFields.size() * 40 + 50;
                    SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
                    SDL_Rect okRect = { panelRect.x + panelRect.w - 190, panelRect.y + panelRect.h - 45, 80, 30 };
                    SDL_Rect cancelRect = { panelRect.x + panelRect.w - 100, panelRect.y + panelRect.h - 45, 80, 30 };
                    if (SDL_PointInRect(&mousePoint, &okRect)) { handle_dialog_ok(); }
                    else if (SDL_PointInRect(&mousePoint, &cancelRect)) { close_dialog(); }
                }
                continue;
            }

            if (currentInteractionMode == InteractionMode::EDITING_COMPONENT_VALUE || currentInteractionMode == InteractionMode::EDITING_LABEL_TEXT) {
                if (e.type == SDL_KEYDOWN) {
                    if (e.key.keysym.sym == SDLK_RETURN || e.key.keysym.sym == SDLK_ESCAPE) {
                        if (e.key.keysym.sym == SDLK_RETURN) {
                            if(editingComponent) editingComponent->value = textInputBuffer;
                            if(editingLabel) editingLabel->text = textInputBuffer;
                        } else if (editingLabel && editingLabel->text.empty()) {
                            labels.pop_back();
                        }
                        editingComponent = nullptr;
                        editingLabel = nullptr;
                        currentInteractionMode = InteractionMode::NONE;
                        SDL_StopTextInput();
                    } else if (e.key.keysym.sym == SDLK_BACKSPACE && !textInputBuffer.empty()) {
                        textInputBuffer.pop_back();
                    } else if (e.key.keysym.sym == SDLK_r && (e.key.keysym.mod & KMOD_CTRL)) {
                        if (editingLabel) {
                            editingLabel->rotation_angle = (editingLabel->rotation_angle + 90) % 360;
                        }
                    }
                } else if (e.type == SDL_TEXTINPUT) {
                    textInputBuffer += e.text.text;
                }
                continue;
            }

            if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_ESCAPE) { currentInteractionMode = InteractionMode::NONE; isDrawingWire = false; }
                else if (e.key.keysym.sym == SDLK_r) {
                    if (currentInteractionMode == InteractionMode::PLACE_COMPONENT) { placementRotation = (placementRotation + 90) % 360; }
                    else if (currentInteractionMode == InteractionMode::DRAGGING_COMPONENT && selectedComponent) {
                        selectedComponent->rotation_angle = (selectedComponent->rotation_angle + 90) % 360;
                        schematicModified = true;
                    }
                }
            }

            if (e.type == SDL_MOUSEBUTTONDOWN) {
                int mouseX, mouseY; SDL_GetMouseState(&mouseX, &mouseY);
                SDL_Point snappedPoint = snap_to_grid(mouseX, mouseY);
                if (e.button.button == SDL_BUTTON_LEFT) {
                    bool valueClicked = false;
                    for (auto& comp : components) {
                        if (SDL_PointInRect(&snappedPoint, &comp.value_rect)) {
                            editingComponent = &comp;
                            textInputBuffer = comp.value;
                            currentInteractionMode = InteractionMode::EDITING_COMPONENT_VALUE;
                            SDL_StartTextInput();
                            valueClicked = true;
                            break;
                        }
                    }
                    if (valueClicked) continue;

                    if (currentInteractionMode == InteractionMode::PLACE_LABEL) {
                        labels.push_back({"", snappedPoint});
                        editingLabel = &labels.back();
                        textInputBuffer = "";
                        currentInteractionMode = InteractionMode::EDITING_LABEL_TEXT;
                        SDL_StartTextInput();
                    } else if (currentInteractionMode == InteractionMode::PLACE_GND_LABEL) {
                        labels.push_back({"GND", snappedPoint, {}, 0});
                        schematicModified = true;
                        currentInteractionMode = InteractionMode::NONE;
                    } else if (currentInteractionMode == InteractionMode::WIRING) {
                        firstInteractionPoint = snappedPoint;
                        isDrawingWire = true;
                    } else if (currentInteractionMode == InteractionMode::PLACE_COMPONENT) {
                        SDL_Point center = snappedPoint;
                        Component newComp;
                        newComp.type = selectedComponentType;
                        newComp.rotation_angle = placementRotation;

                        int half_len = COMPONENT_DEFAULT_LENGTH / 2;
                        if (placementRotation == 0) { newComp.node1 = {center.x - half_len, center.y}; newComp.node2 = {center.x + half_len, center.y}; }
                        else if (placementRotation == 90) { newComp.node1 = {center.x, center.y - half_len}; newComp.node2 = {center.x, center.y + half_len}; }
                        else if (placementRotation == 180) { newComp.node1 = {center.x + half_len, center.y}; newComp.node2 = {center.x - half_len, center.y}; }
                        else { newComp.node1 = {center.x, center.y + half_len}; newComp.node2 = {center.x, center.y - half_len}; }

                        newComp.id = generate_component_name(newComp.type);
                        get_component_defaults(newComp.type, newComp.value);
                        components.push_back(newComp);
                        schematicModified = true;
                    } else if (currentInteractionMode == InteractionMode::DELETE_ITEM) {
                        const double DELETE_THRESHOLD_SQ = 10 * 10;
                        SDL_Point clickPoint = {mouseX, mouseY};

                        double min_dist_sq = DELETE_THRESHOLD_SQ;
                        int wire_to_delete = -1;
                        int component_to_delete = -1;
                        int label_to_delete = -1;

                        for (int i = 0; i < labels.size(); ++i) {
                            if (SDL_PointInRect(&clickPoint, &labels[i].text_rect)) {
                                label_to_delete = i;
                                break;
                            }
                        }

                        if (label_to_delete == -1) {
                            for (int i = 0; i < wires.size(); ++i) {
                                double dist_sq = dist_to_segment_sq(clickPoint, wires[i].start, wires[i].end);
                                if (dist_sq < min_dist_sq) {
                                    min_dist_sq = dist_sq;
                                    wire_to_delete = i;
                                    component_to_delete = -1;
                                }
                            }
                            for (int i = 0; i < components.size(); ++i) {
                                double dist_sq = dist_to_segment_sq(clickPoint, components[i].node1, components[i].node2);
                                if (dist_sq < min_dist_sq) {
                                    min_dist_sq = dist_sq;
                                    component_to_delete = i;
                                    wire_to_delete = -1;
                                }
                            }
                        }

                        if (label_to_delete != -1) {
                            labels.erase(labels.begin() + label_to_delete);
                            schematicModified = true;
                        } else if (wire_to_delete != -1) {
                            wires.erase(wires.begin() + wire_to_delete);
                            schematicModified = true;
                        } else if (component_to_delete != -1) {
                            components.erase(components.begin() + component_to_delete);
                            schematicModified = true;
                        }
                    }
                } else if (e.button.button == SDL_BUTTON_RIGHT) {
                    currentInteractionMode = InteractionMode::NONE;
                    isDrawingWire = false;
                }
            } else if (e.type == SDL_MOUSEBUTTONUP) {
                if (e.button.button == SDL_BUTTON_LEFT) {
                    if (currentInteractionMode == InteractionMode::WIRING && isDrawingWire) {
                        int mouseX, mouseY; SDL_GetMouseState(&mouseX, &myY);
                        SDL_Point snappedPoint = snap_to_grid(mouseX, mouseY);
                        if (snappedPoint.x != firstInteractionPoint.x || snappedPoint.y != firstInteractionPoint.y) {
                            wires.push_back({firstInteractionPoint, snappedPoint});
                            schematicModified = true;
                        }
                        isDrawingWire = false;
                    }
                }
            }

            if (fileMenuButtons[5].m_is_hovered && e.type == SDL_MOUSEBUTTONDOWN) quit = true;
            if (isFileMenuOpen) for (auto& b : fileMenuButtons) b.handle_event(&e);
            if (isComponentsMenuOpen) for (auto& i : componentMenuItems) i.button.handle_event(&e);
            if (isSimulateMenuOpen) for (auto& b : simulateMenuButtons) b.handle_event(&e);
            for (auto& b : topBarButtons) b.handle_event(&e);
        }

        SDL_SetRenderDrawColor(renderer, 0x22, 0x22, 0x22, 0xFF);
        SDL_RenderClear(renderer);

        draw_grid(renderer);
        draw_schematic_elements(renderer, valueFont, componentMenuItems);

        if (currentInteractionMode == InteractionMode::WIRING && isDrawingWire) {
            int mx, my; SDL_GetMouseState(&mx, &my);
            SDL_Point currentSnappedPoint = snap_to_grid(mx, my);
            thickLineRGBA(renderer, firstInteractionPoint.x, firstInteractionPoint.y, currentSnappedPoint.x, currentSnappedPoint.y, 3, 0xFF, 0xFF, 0xFF, 0xFF);
        }

        if (currentInteractionMode == InteractionMode::PLACE_COMPONENT) {
            int mx, my; SDL_GetMouseState(&mx, &my);
            SDL_Point center = snap_to_grid(mx, my);
            SDL_Texture* icon = nullptr;
            for(const auto& item : componentMenuItems) if(item.type == selectedComponentType) icon = item.iconTexture;
            if (icon) {
                SDL_SetTextureAlphaMod(icon, 150);
                int icon_size = COMPONENT_DEFAULT_LENGTH;
                int icon_half_size = icon_size / 2;
                SDL_Rect dest = {center.x - icon_half_size, center.y - icon_half_size, icon_size, icon_size};
                SDL_Point rotationCenter = { icon_half_size, icon_half_size };
                SDL_RenderCopyEx(renderer, icon, nullptr, &dest, placementRotation, &rotationCenter, SDL_FLIP_NONE);
                SDL_SetTextureAlphaMod(icon, 255);
            }
        }

        SDL_Rect topBarRect = {0, 0, SCREEN_WIDTH, TOP_BAR_HEIGHT};
        SDL_SetRenderDrawColor(renderer, 0x33, 0x33, 0x33, 0xFF);
        SDL_RenderFillRect(renderer, &topBarRect);
        for (auto& b : topBarButtons) b.render(renderer, uiFont);
        if (isFileMenuOpen) for (auto& b : fileMenuButtons) b.render(renderer, uiFont);
        if (isSimulateMenuOpen) for (auto& b : simulateMenuButtons) b.render(renderer, uiFont);
        if (isComponentsMenuOpen) {
            SDL_Rect menuPanel = { 205, 35, 4*(COMP_ITEM_W + COMP_ITEM_PAD) + 5, 2*(COMP_ITEM_H + COMP_ITEM_PAD) + 5};
            SDL_SetRenderDrawColor(renderer, 0x3A, 0x3A, 0x3A, 0xFF);
            SDL_RenderFillRect(renderer, &menuPanel);
            for (auto& item : componentMenuItems) {
                item.button.render(renderer, uiFont);
                if (item.iconTexture) {
                    SDL_Rect iconRect = { item.button.m_position.x + 45, item.button.m_position.y + 10, 60, 60 };
                    SDL_RenderCopy(renderer, item.iconTexture, nullptr, &iconRect);
                }
                SDL_Color textColor = { 0xFF, 0xFF, 0xFF, 0xFF };
                SDL_Surface* textSurface = TTF_RenderText_Solid(uiFont, item.name.c_str(), textColor);
                if (textSurface) {
                    SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                    SDL_Rect textRect = { item.button.m_position.x + (item.button.m_position.w - textSurface->w) / 2, item.button.m_position.y + 75, textSurface->w, textSurface->h };
                    SDL_RenderCopy(renderer, textTexture, nullptr, &textRect);
                    SDL_FreeSurface(textSurface);
                    SDL_DestroyTexture(textTexture);
                }
            }
        }

        if (currentInteractionMode == InteractionMode::DIALOG_ACTIVE) {
            render_dialog(renderer, uiFont);
        }

        SDL_RenderPresent(renderer);
    }

    SDL_StopTextInput();
    for(auto& item : componentMenuItems) { if(item.iconTexture) SDL_DestroyTexture(item.iconTexture); }
    TTF_CloseFont(uiFont);
    TTF_CloseFont(valueFont);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    IMG_Quit();
    TTF_Quit();
    SDL_Quit();
    return 0;
}
