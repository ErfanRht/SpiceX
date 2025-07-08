#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>
#ifdef _WIN32
#include <SDL2/SDL2_gfx.h>
#else
#include <SDL2/SDL2_gfxPrimitives.h>
#endif
#include <bits/stdc++.h>

using namespace std;

const int SCREEN_WIDTH = 1280;
const int SCREEN_HEIGHT = 720;
const int TOP_BAR_HEIGHT = 40;
const int GRID_SPACING = 20;

bool isFileMenuOpen = false;
bool isComponentsMenuOpen = false;
bool isSimulateMenuOpen = false;


class Button {
public:
    Button(string text, int x, int y, int w, int h, function<void()> onClick) {
        m_text = text;
        m_position = {x, y, w, h};
        m_onClick = onClick;
        m_is_hovered = false;
    }

    SDL_Rect m_position;

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
                int text_w = textSurface->w;
                int text_h = textSurface->h;
                SDL_FreeSurface(textSurface);
                SDL_Rect renderQuad = { m_position.x + (m_position.w - text_w) / 2, m_position.y + (m_position.h - text_h) / 2, text_w, text_h };
                SDL_RenderCopy(renderer, textTexture, nullptr, &renderQuad);
                SDL_DestroyTexture(textTexture);
            }
        }
    }

private:
    string m_text;
    bool m_is_hovered;
    function<void()> m_onClick;
};

struct ComponentMenuItem {
    string name;
    SDL_Texture* iconTexture = nullptr;
    Button button;
};


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
    if (originalSurface == nullptr) {
        cerr << "Failed to load image " << path << "! SDL_image Error: " << IMG_GetError() << endl;
        return nullptr;
    }

    SDL_LockSurface(originalSurface);
    uint32_t* pixels = (uint32_t*)originalSurface->pixels;
    int pixelCount = originalSurface->w * originalSurface->h;
    SDL_PixelFormat* format = originalSurface->format;
    uint32_t white = SDL_MapRGB(format, 255, 255, 255);

    for (int i = 0; i < pixelCount; ++i) {
        uint8_t r, g, b, a;
        SDL_GetRGBA(pixels[i], format, &r, &g, &b, &a);
        if (a > 50) {
            pixels[i] = white;
        }
    }

    SDL_UnlockSurface(originalSurface);
    SDL_Texture* finalTexture = SDL_CreateTextureFromSurface(renderer, originalSurface);
    SDL_FreeSurface(originalSurface);
    SDL_SetTextureBlendMode(finalTexture, SDL_BLENDMODE_BLEND);

    return finalTexture;
}


int main(int argc, char* args[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0 || TTF_Init() == -1) { cerr << "Core SDL init failed." << endl; return -1; }
    int imgFlags = IMG_INIT_PNG;
    if (!(IMG_Init(imgFlags) & imgFlags)) { cerr << "SDL_image init failed." << endl; return -1; }

    SDL_Window* window = SDL_CreateWindow("SUTSpice | Phase 2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    TTF_Font* font = TTF_OpenFont(R"(C:\Windows\Fonts\arial.ttf)", 16);
    if (!window || !renderer || !font) { cerr << "Window, renderer, or font creation failed." << endl; return -1; }

    const string ASSET_PATH = "C:/Users/Erfan/Dev/Cpp/sutSpice2/assets/";

    vector<Button> topBarButtons;
    topBarButtons.emplace_back("File", 10, 5, 80, 30, [](){ isFileMenuOpen = !isFileMenuOpen; isComponentsMenuOpen = false; isSimulateMenuOpen = false; });
    topBarButtons.emplace_back("Simulate", 100, 5, 100, 30, [](){ isSimulateMenuOpen = !isSimulateMenuOpen; isFileMenuOpen = false; isComponentsMenuOpen = false; });
    topBarButtons.emplace_back("Add Component", 210, 5, 150, 30, [](){ isComponentsMenuOpen = !isComponentsMenuOpen; isFileMenuOpen = false; isSimulateMenuOpen = false; });

    vector<Button> fileMenuButtons;
    fileMenuButtons.emplace_back("New Schematic", 10, 40, 180, 30, [](){ cout << "Action: New" << endl; isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Open Schematic", 10, 70, 180, 30, [](){ cout << "Action: Open" << endl; isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Save", 10, 100, 180, 30, [](){ cout << "Action: Save" << endl; isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Save As...", 10, 130, 180, 30, [](){ cout << "Action: Save As" << endl; isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Clear", 10, 160, 180, 30, [](){ cout << "Action: Clear" << endl; isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Exit", 10, 190, 180, 30, [](){ isFileMenuOpen = false; });

    vector<Button> simulateMenuButtons;
    simulateMenuButtons.emplace_back("DC Sweep Analysis", 100, 40, 180, 30, [](){ cout << "Action: DC Sweep" << endl; isSimulateMenuOpen = false; });
    simulateMenuButtons.emplace_back("Transient Analysis", 100, 70, 180, 30, [](){ cout << "Action: Transient" << endl; isSimulateMenuOpen = false; });

    vector<ComponentMenuItem> componentMenuItems;
    const int COMPONENT_MENU_START_X = 210;
    const int COMPONENT_MENU_START_Y = 40;
    const int COMPONENT_ITEM_WIDTH = 150;
    const int COMPONENT_ITEM_HEIGHT = 100;
    const int COMPONENT_ITEM_PADDING = 10;

    componentMenuItems.push_back({"Resistor", loadAndProcessTexture(ASSET_PATH + "resistor.png", renderer), Button("", COMPONENT_MENU_START_X, COMPONENT_MENU_START_Y, COMPONENT_ITEM_WIDTH, COMPONENT_ITEM_HEIGHT, [](){ cout << "Add Resistor" << endl; isComponentsMenuOpen = false; })});
    componentMenuItems.push_back({"Capacitor", loadAndProcessTexture(ASSET_PATH + "capacitor.png", renderer), Button("", COMPONENT_MENU_START_X + 1*(COMPONENT_ITEM_WIDTH + COMPONENT_ITEM_PADDING), COMPONENT_MENU_START_Y, COMPONENT_ITEM_WIDTH, COMPONENT_ITEM_HEIGHT, [](){ cout << "Add Capacitor" << endl; isComponentsMenuOpen = false; })});
    componentMenuItems.push_back({"Inductor", loadAndProcessTexture(ASSET_PATH + "inductor.png", renderer), Button("", COMPONENT_MENU_START_X + 2*(COMPONENT_ITEM_WIDTH + COMPONENT_ITEM_PADDING), COMPONENT_MENU_START_Y, COMPONENT_ITEM_WIDTH, COMPONENT_ITEM_HEIGHT, [](){ cout << "Add Inductor" << endl; isComponentsMenuOpen = false; })});
    componentMenuItems.push_back({"Diode", loadAndProcessTexture(ASSET_PATH + "diode.png", renderer), Button("", COMPONENT_MENU_START_X + 3*(COMPONENT_ITEM_WIDTH + COMPONENT_ITEM_PADDING), COMPONENT_MENU_START_Y, COMPONENT_ITEM_WIDTH, COMPONENT_ITEM_HEIGHT, [](){ cout << "Add Diode" << endl; isComponentsMenuOpen = false; })});
    componentMenuItems.push_back({"DC Voltage Source", loadAndProcessTexture(ASSET_PATH + "voltage_source.png", renderer), Button("", COMPONENT_MENU_START_X, COMPONENT_MENU_START_Y + 1*(COMPONENT_ITEM_HEIGHT + COMPONENT_ITEM_PADDING), COMPONENT_ITEM_WIDTH, COMPONENT_ITEM_HEIGHT, [](){ cout << "Add DC Voltage Source" << endl; isComponentsMenuOpen = false; })});
    componentMenuItems.push_back({"DC Current Source", loadAndProcessTexture(ASSET_PATH + "current_source.png", renderer), Button("", COMPONENT_MENU_START_X + 1*(COMPONENT_ITEM_WIDTH + COMPONENT_ITEM_PADDING), COMPONENT_MENU_START_Y + 1*(COMPONENT_ITEM_HEIGHT + COMPONENT_ITEM_PADDING), COMPONENT_ITEM_WIDTH, COMPONENT_ITEM_HEIGHT, [](){ cout << "Add DC Current Source" << endl; isComponentsMenuOpen = false; })});
    componentMenuItems.push_back({"AC Voltage Source", loadAndProcessTexture(ASSET_PATH + "ac_voltage_source.png", renderer), Button("", COMPONENT_MENU_START_X + 2*(COMPONENT_ITEM_WIDTH + COMPONENT_ITEM_PADDING), COMPONENT_MENU_START_Y + 1*(COMPONENT_ITEM_HEIGHT + COMPONENT_ITEM_PADDING), COMPONENT_ITEM_WIDTH, COMPONENT_ITEM_HEIGHT, [](){ cout << "Add AC Voltage Source" << endl; isComponentsMenuOpen = false; })});

    bool quit = false;
    SDL_Event e;
    while (!quit) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) { quit = true; }
            if (e.type == SDL_MOUSEBUTTONDOWN) {
                int x, y;
                SDL_GetMouseState(&x, &y);
                if(isFileMenuOpen && x > 10 && x < 190 && y > 190 && y < 220) quit = true;
            }

            if (isFileMenuOpen) { for (auto& button : fileMenuButtons) button.handle_event(&e); }
            if (isComponentsMenuOpen) { for (auto& item : componentMenuItems) item.button.handle_event(&e); }
            if (isSimulateMenuOpen) { for (auto& button : simulateMenuButtons) button.handle_event(&e); }

            for (auto& button : topBarButtons) button.handle_event(&e);
        }

        SDL_SetRenderDrawColor(renderer, 0x22, 0x22, 0x22, 0xFF);
        SDL_RenderClear(renderer);
        draw_grid(renderer);

        SDL_Rect top_bar_rect = { 0, 0, SCREEN_WIDTH, TOP_BAR_HEIGHT };
        SDL_SetRenderDrawColor(renderer, 0x33, 0x33, 0x33, 0xFF);
        SDL_RenderFillRect(renderer, &top_bar_rect);
        for (auto& button : topBarButtons) button.render(renderer, font);

        if (isFileMenuOpen) {
            for (auto& button : fileMenuButtons) button.render(renderer, font);
        }

        if (isComponentsMenuOpen) {
            SDL_Rect menuPanel = { COMPONENT_MENU_START_X - 5, COMPONENT_MENU_START_Y - 5, 4 * (COMPONENT_ITEM_WIDTH + COMPONENT_ITEM_PADDING) + 5, 2 * (COMPONENT_ITEM_HEIGHT + COMPONENT_ITEM_PADDING) + 5};
            SDL_SetRenderDrawColor(renderer, 0x3A, 0x3A, 0x3A, 0xFF);
            SDL_RenderFillRect(renderer, &menuPanel);

            for (auto& item : componentMenuItems) {
                item.button.render(renderer, font);
                if (item.iconTexture) {
                    SDL_Rect iconRect = { item.button.m_position.x + 45, item.button.m_position.y + 10, 60, 60 };
                    SDL_RenderCopy(renderer, item.iconTexture, nullptr, &iconRect);
                }
                SDL_Color textColor = { 0xFF, 0xFF, 0xFF, 0xFF };
                SDL_Surface* textSurface = TTF_RenderText_Solid(font, item.name.c_str(), textColor);
                if (textSurface) {
                    SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                    int text_w = textSurface->w;
                    int text_h = textSurface->h;
                    SDL_FreeSurface(textSurface);
                    SDL_Rect textRect = { item.button.m_position.x + (item.button.m_position.w - text_w) / 2, item.button.m_position.y + 75, text_w, text_h };
                    SDL_RenderCopy(renderer, textTexture, nullptr, &textRect);
                    SDL_DestroyTexture(textTexture);
                }
            }
        }

        if (isSimulateMenuOpen) {
            for (auto& button : simulateMenuButtons) button.render(renderer, font);
        }

        SDL_RenderPresent(renderer);
    }

    for(auto& item : componentMenuItems) { if(item.iconTexture) SDL_DestroyTexture(item.iconTexture); }
    TTF_CloseFont(font);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    IMG_Quit();
    TTF_Quit();
    SDL_Quit();

    return 0;
}
