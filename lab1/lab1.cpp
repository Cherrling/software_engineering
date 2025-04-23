#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <limits>  // 添加到其他头文件后面
#include <ctime>  // 添加到其他头文件后面
#include <queue>  // 添加到其他头文件后面
#include <climits>  // 添加到其他头文件后面，用于INT_MAX
// 定义有向图的数据结构
struct Edge {
    std::string to;
    int weight;
    
    Edge(const std::string& to, int weight) : to(to), weight(weight) {}
};
// 定义用于Dijkstra算法的路径结构
struct PathInfo {
    int distance;
    std::vector<std::string> path;
    
    PathInfo() : distance(INT_MAX) {} // 初始化为无穷大
    
    PathInfo(int d, const std::vector<std::string>& p) : distance(d), path(p) {}
};
// 转换为小写
std::string toLower(const std::string& word) {
    std::string result = word;
    std::transform(result.begin(), result.end(), result.begin(), 
                [](unsigned char c) { return std::tolower(c); });
    return result;
}
class DirectedGraph {
private:
    std::map<std::string, std::vector<Edge>> adjacencyList;

public:
    // 添加或更新边
    void addEdge(const std::string& from, const std::string& to) {
        // 查找是否已存在该边
        for (auto& edge : adjacencyList[from]) {
            if (edge.to == to) {
                // 边已存在，增加权重
                edge.weight++;
                return;
            }
        }
        // 边不存在，创建新边
        adjacencyList[from].push_back(Edge(to, 1));
    }

    // 优化CLI界面展示
    void printGraph() {
        std::cout << "+" << std::string(60, '-') << "+" << std::endl;
        std::cout << "|" << std::setw(30) << "有向图结构" << std::setw(30) << "|" << std::endl;
        std::cout << "+" << std::string(60, '-') << "+" << std::endl;
        
        for (const auto& node : adjacencyList) {
            std::cout << "| " << std::setw(10) << std::left << node.first << " | ";
            
            if (node.second.empty()) {
                std::cout << std::setw(46) << std::left << "无出边" << "|" << std::endl;
            } else {
                bool firstEdge = true;
                for (const auto& edge : node.second) {
                    if (!firstEdge) {   
                        std::cout << "| " << std::setw(10) << "" << " | ";
                    }
                    
                    std::string edgeStr = "-> " + edge.to + " (权重=" + std::to_string(edge.weight) + ")";
                    std::cout << std::setw(46) << std::left << edgeStr << "|" << std::endl;
                    firstEdge = false;
                }
            }
            std::cout << "+" << std::string(60, '-') << "+" << std::endl;
        }
        
        std::cout << "+" << std::string(60, '-') << "+" << std::endl;
    }
    
    // 以DOT格式导出图形
    bool exportToDOT(const std::string& filename) {
        std::ofstream dotFile(filename);
        if (!dotFile.is_open()) {
            std::cerr << "无法创建DOT文件: " << filename << std::endl;
            return false;
        }
        
        dotFile << "digraph TextGraph {" << std::endl;
        dotFile << "    rankdir=LR;" << std::endl;
        dotFile << "    node [shape=box, style=filled, fillcolor=lightblue];" << std::endl;
        dotFile << "    edge [color=darkblue];" << std::endl;
        
        // 添加所有节点
        for (const auto& node : adjacencyList) {
            dotFile << "    \"" << node.first << "\" [label=\"" << node.first << "\"];" << std::endl;
        }
        
        // 添加所有边
        for (const auto& node : adjacencyList) {
            for (const auto& edge : node.second) {
                dotFile << "    \"" << node.first << "\" -> \"" << edge.to 
                       << "\" [label=\"" << edge.weight << "\", penwidth=" 
                       << std::min(1 + edge.weight * 0.5, 5.0) << "];" << std::endl;
            }
        }
        
        dotFile << "}" << std::endl;
        dotFile.close();
        
        return true;
    }
    
    // 使用Graphviz生成图像
    bool generateImage(const std::string& dotFilename, const std::string& outputFilename) {
        std::string command = "dot -Tpng " + dotFilename + " -o " + outputFilename;
        int result = system(command.c_str());
        
        if (result != 0) {
            std::cerr << "生成图像失败。请确保已安装Graphviz。" << std::endl;
            std::cerr << "可以通过运行 'sudo apt-get install graphviz' 安装Graphviz。" << std::endl;
            return false;
        }
        
        std::cout << "图像已成功生成: " << outputFilename << std::endl;
        return true;
    }

    
    // 查找桥接词
    std::vector<std::string> findBridgeWords(const std::string& word1, const std::string& word2) {
        std::vector<std::string> bridgeWords;
        
        // 检查两个词是否都在图中
        if (adjacencyList.find(word1) == adjacencyList.end() || 
            adjacencyList.find(word2) == adjacencyList.end()) {
            return bridgeWords; // 空列表表示word1或word2不在图中
        }
        
        // 查找从word1出发的所有边
        for (const auto& edge : adjacencyList[word1]) {
            std::string potentialBridge = edge.to;
            
            // 检查这个潜在的桥接词是否有一条边指向word2
            if (adjacencyList.find(potentialBridge) != adjacencyList.end()) {
                for (const auto& nextEdge : adjacencyList[potentialBridge]) {
                    if (nextEdge.to == word2) {
                        // 找到桥接词
                        bridgeWords.push_back(potentialBridge);
                        break;
                    }
                }
            }
        }
        
        return bridgeWords;
    }
    
    // 检查单词是否在图中
    bool containsWord(const std::string& word) {
        return adjacencyList.find(word) != adjacencyList.end();
    }

    
    // 生成新文本，插入桥接词
    std::string generateNewText(const std::string& inputText) {
        std::vector<std::string> words;
        std::string currentWord;
        
        // 分词
        for (char ch : inputText) {
            if (std::isalpha(ch)) {
                currentWord += ch;
            } else if (!currentWord.empty()) {
                words.push_back(toLower(currentWord));
                currentWord = "";
            }
        }
        
        // 处理最后一个单词
        if (!currentWord.empty()) {
            words.push_back(toLower(currentWord));
        }
        
        // 如果少于2个单词，直接返回原文本
        if (words.size() < 2) {
            return inputText;
        }
        
        // 构建新文本
        std::string newText = words[0];
        
        for (size_t i = 0; i < words.size() - 1; ++i) {
            const std::string& word1 = words[i];
            const std::string& word2 = words[i + 1];
            
            // 查找桥接词
            std::vector<std::string> bridges = findBridgeWords(word1, word2);
            
            if (!bridges.empty()) {
                // 随机选择一个桥接词
                srand(static_cast<unsigned int>(time(nullptr)) + i); // 为每对单词使用不同的种子
                std::string selectedBridge = bridges[rand() % bridges.size()];
                
                newText += " " + selectedBridge;
            }
            
            newText += " " + word2;
        }
        
        return newText;
    }



    // 计算两个单词之间的最短路径
    PathInfo findShortestPath(const std::string& start, const std::string& end) {
        // 检查开始和结束单词是否在图中
        if (adjacencyList.find(start) == adjacencyList.end() || 
            adjacencyList.find(end) == adjacencyList.end()) {
            return PathInfo(); // 返回无效路径
        }
        
        // 初始化距离表和前驱表
        std::map<std::string, int> distance;
        std::map<std::string, std::string> previous;
        std::map<std::string, bool> visited;
        
        for (const auto& node : adjacencyList) {
            distance[node.first] = INT_MAX;
            visited[node.first] = false;
        }
        
        // 设置起始点距离为0
        distance[start] = 0;
        
        // Dijkstra算法
        while (true) {
            // 找到未访问节点中距离最小的
            std::string current = "";
            int minDist = INT_MAX;
            
            for (const auto& node : adjacencyList) {
                if (!visited[node.first] && distance[node.first] < minDist) {
                    minDist = distance[node.first];
                    current = node.first;
                }
            }
            
            // 如果没有找到有效节点或者已到达终点，结束循环
            if (current.empty() || current == end) break;
            
            // 标记为已访问
            visited[current] = true;
            
            // 更新相邻节点的距离
            for (const auto& edge : adjacencyList[current]) {
                int newDist = distance[current] + edge.weight;
                
                if (newDist < distance[edge.to]) {
                    distance[edge.to] = newDist;
                    previous[edge.to] = current;
                }
            }
        }
        
        // 如果终点不可达
        if (distance[end] == INT_MAX) {
            return PathInfo(); // 返回无效路径
        }
        
        // 重建路径
        std::vector<std::string> path;
        for (std::string at = end; at != start; at = previous[at]) {
            path.push_back(at);
        }
        path.push_back(start);
        
        // 反转路径，使其从起点开始
        std::reverse(path.begin(), path.end());
        
        return PathInfo(distance[end], path);
    }
    
    // 生成带高亮路径的DOT文件
    bool exportToDOTWithPath(const std::string& filename, const std::vector<std::string>& path) {
        std::ofstream dotFile(filename);
        if (!dotFile.is_open()) {
            std::cerr << "无法创建DOT文件: " << filename << std::endl;
            return false;
        }
        
        // 创建边的路径映射，用于快速检查边是否在路径中
        std::map<std::string, std::map<std::string, bool>> pathEdges;
        
        for (size_t i = 0; i < path.size() - 1; ++i) {
            pathEdges[path[i]][path[i + 1]] = true;
        }
        
        dotFile << "digraph TextGraph {" << std::endl;
        dotFile << "    rankdir=LR;" << std::endl;
        dotFile << "    node [shape=box];" << std::endl;
        
        // 添加所有节点，高亮路径中的节点
        for (const auto& node : adjacencyList) {
            std::string nodeStyle = "style=filled, ";
            
            // 检查节点是否在路径中
            bool inPath = (std::find(path.begin(), path.end(), node.first) != path.end());
            
            if (inPath) {
                // 路径上的节点使用不同的颜色
                nodeStyle += "fillcolor=orange";
            } else {
                nodeStyle += "fillcolor=lightblue";
            }
            
            // 如果是起点或终点，使用特别的标记
            if (!path.empty() && (node.first == path.front() || node.first == path.back())) {
                nodeStyle += ", penwidth=2";
            }
            
            dotFile << "    \"" << node.first << "\" [" << nodeStyle << "];" << std::endl;
        }
        
        // 添加所有边，高亮路径中的边
        for (const auto& node : adjacencyList) {
            for (const auto& edge : node.second) {
                std::string edgeStyle = "";
                
                // 检查边是否在路径中
                bool inPath = (pathEdges.find(node.first) != pathEdges.end() && 
                             pathEdges[node.first].find(edge.to) != pathEdges[node.first].end());
                
                if (inPath) {
                    // 路径上的边使用不同的颜色和粗细
                    edgeStyle = "color=red, penwidth=2";
                } else {
                    edgeStyle = "color=darkblue";
                }
                
                dotFile << "    \"" << node.first << "\" -> \"" << edge.to 
                       << "\" [label=\"" << edge.weight << "\", " << edgeStyle << "];" << std::endl;
            }
        }
        
        dotFile << "}" << std::endl;
        dotFile.close();
        
        return true;
    }
    
    // 计算一个单词到所有其他单词的最短路径
    std::map<std::string, PathInfo> findAllShortestPaths(const std::string& start) {
        std::map<std::string, PathInfo> allPaths;
        
        // 检查开始单词是否在图中
        if (adjacencyList.find(start) == adjacencyList.end()) {
            return allPaths; // 返回空映射
        }
        
        // 初始化距离表和前驱表
        std::map<std::string, int> distance;
        std::map<std::string, std::string> previous;
        std::map<std::string, bool> visited;
        
        for (const auto& node : adjacencyList) {
            distance[node.first] = INT_MAX;
            visited[node.first] = false;
        }
        
        // 设置起始点距离为0
        distance[start] = 0;
        
        // Dijkstra算法
        while (true) {
            // 找到未访问节点中距离最小的
            std::string current = "";
            int minDist = INT_MAX;
            
            for (const auto& node : adjacencyList) {
                if (!visited[node.first] && distance[node.first] < minDist) {
                    minDist = distance[node.first];
                    current = node.first;
                }
            }
            
            // 如果没有找到有效节点，结束循环
            if (current.empty() || minDist == INT_MAX) break;
            
            // 标记为已访问
            visited[current] = true;
            
            // 更新相邻节点的距离
            for (const auto& edge : adjacencyList[current]) {
                int newDist = distance[current] + edge.weight;
                
                if (newDist < distance[edge.to]) {
                    distance[edge.to] = newDist;
                    previous[edge.to] = current;
                }
            }
        }
        
        // 为每个可达节点构建路径信息
        for (const auto& node : adjacencyList) {
            if (node.first != start && distance[node.first] != INT_MAX) {
                std::vector<std::string> path;
                for (std::string at = node.first; at != start; at = previous[at]) {
                    path.push_back(at);
                }
                path.push_back(start);
                
                // 反转路径，使其从起点开始
                std::reverse(path.begin(), path.end());
                
                allPaths[node.first] = PathInfo(distance[node.first], path);
            }
        }
        
        return allPaths;
    }
    // 在DirectedGraph类中添加PageRank算法实现

    // 计算图中所有节点的PageRank值
    std::map<std::string, double> calculatePageRank(double dampingFactor = 0.85, int iterations = 100, double tolerance = 1e-6) {
        // 初始化PageRank值
        std::map<std::string, double> pageRank;
        std::map<std::string, double> newPageRank;
        double initialRank = 1.0 / adjacencyList.size(); // 均匀初始化
        
        for (const auto& node : adjacencyList) {
            pageRank[node.first] = initialRank;
        }
        
        // 迭代计算PageRank值
        for (int i = 0; i < iterations; ++i) {
            // 随机游走的概率部分
            double randomJumpProbability = (1.0 - dampingFactor) / adjacencyList.size();
            
            // 初始化新的PageRank值为随机游走概率
            for (const auto& node : adjacencyList) {
                newPageRank[node.first] = randomJumpProbability;
            }
            
            // 计算链接贡献
            for (const auto& node : adjacencyList) {
                std::string sourceNode = node.first;
                const std::vector<Edge>& outLinks = node.second;
                
                // 如果节点没有出边，则平均分配给所有节点
                if (outLinks.empty()) {
                    double distributedRank = dampingFactor * pageRank[sourceNode] / adjacencyList.size();
                    for (const auto& targetNode : adjacencyList) {
                        newPageRank[targetNode.first] += distributedRank;
                    }
                } else {
                    // 计算当前节点出边的总权重
                    int totalWeight = 0;
                    for (const auto& edge : outLinks) {
                        totalWeight += edge.weight;
                    }
                    
                    // 根据权重分配PageRank值
                    for (const auto& edge : outLinks) {
                        double contribution = dampingFactor * pageRank[sourceNode] * edge.weight / totalWeight;
                        newPageRank[edge.to] += contribution;
                    }
                }
            }
            
            // 检查收敛性
            double diff = 0.0;
            for (const auto& node : adjacencyList) {
                diff += std::abs(newPageRank[node.first] - pageRank[node.first]);
            }
            
            // 更新PageRank值
            pageRank = newPageRank;
            
            // 如果变化小于阈值，则认为已收敛
            if (diff < tolerance) {
                std::cout << "PageRank收敛于迭代" << i + 1 << "（总迭代次数：" << iterations << "）" << std::endl;
                break;
            }
        }
        
        return pageRank;
    }

    // 使用TF-IDF方式改进初始PageRank值
    std::map<std::string, double> calculatePageRankWithTFIDF(double dampingFactor = 0.85, int iterations = 100, double tolerance = 1e-6) {
        // 计算每个节点的出度和入度
        std::map<std::string, int> outDegree;
        std::map<std::string, int> inDegree;
        
        for (const auto& node : adjacencyList) {
            outDegree[node.first] = node.second.size();
            
            for (const auto& edge : node.second) {
                inDegree[edge.to]++;
            }
        }
        
        // 初始化PageRank值，使用归一化的入度作为初始值
        // 这模拟了TF-IDF的思想：被引用越多的节点越重要
        std::map<std::string, double> pageRank;
        std::map<std::string, double> newPageRank;
        double totalInDegree = 0;
        
        for (const auto& degree : inDegree) {
            totalInDegree += degree.second;
        }
        
        // 如果没有入度，使用均匀分布
        if (totalInDegree == 0) {
            double initialRank = 1.0 / adjacencyList.size();
            for (const auto& node : adjacencyList) {
                pageRank[node.first] = initialRank;
            }
        } else {
            // 使用归一化的入度作为初始PR值
            double minPR = 0.5 / adjacencyList.size(); // 设置最小PR值，避免完全没有入边的节点PR为0
            
            for (const auto& node : adjacencyList) {
                pageRank[node.first] = minPR;
                if (inDegree.find(node.first) != inDegree.end()) {
                    pageRank[node.first] += 0.5 * inDegree[node.first] / totalInDegree;
                }
            }
            
            // 归一化初始PR值总和为1
            double sumPR = 0;
            for (const auto& pr : pageRank) {
                sumPR += pr.second;
            }
            
            for (auto& pr : pageRank) {
                pr.second /= sumPR;
            }
        }
        
        // 迭代计算PageRank值（与普通PageRank相同）
        for (int i = 0; i < iterations; ++i) {
            // 随机游走的概率部分
            double randomJumpProbability = (1.0 - dampingFactor) / adjacencyList.size();
            
            // 初始化新的PageRank值为随机游走概率
            for (const auto& node : adjacencyList) {
                newPageRank[node.first] = randomJumpProbability;
            }
            
            // 计算链接贡献
            for (const auto& node : adjacencyList) {
                std::string sourceNode = node.first;
                const std::vector<Edge>& outLinks = node.second;
                
                // 如果节点没有出边，则平均分配给所有节点
                if (outLinks.empty()) {
                    double distributedRank = dampingFactor * pageRank[sourceNode] / adjacencyList.size();
                    for (const auto& targetNode : adjacencyList) {
                        newPageRank[targetNode.first] += distributedRank;
                    }
                } else {
                    // 计算当前节点出边的总权重
                    int totalWeight = 0;
                    for (const auto& edge : outLinks) {
                        totalWeight += edge.weight;
                    }
                    
                    // 根据权重分配PageRank值
                    for (const auto& edge : outLinks) {
                        double contribution = dampingFactor * pageRank[sourceNode] * edge.weight / totalWeight;
                        newPageRank[edge.to] += contribution;
                    }
                }
            }
            
            // 检查收敛性
            double diff = 0.0;
            for (const auto& node : adjacencyList) {
                diff += std::abs(newPageRank[node.first] - pageRank[node.first]);
            }
            
            // 更新PageRank值
            pageRank = newPageRank;
            
            // 如果变化小于阈值，则认为已收敛
            if (diff < tolerance) {
                std::cout << "PageRank收敛于迭代" << i + 1 << "（总迭代次数：" << iterations << "）" << std::endl;
                break;
            }
        }
        
        return pageRank;
    }
    // 在DirectedGraph类中添加随机游走功能

    // 随机游走
    std::vector<std::string> randomWalk() {
        // 确保图非空
        if (adjacencyList.empty()) {
            return {};
        }
        // 清空缓冲区
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        
        // 随机选择起始节点
        std::srand(std::time(nullptr));
        int startNodeIndex = std::rand() % adjacencyList.size();
        
        auto it = adjacencyList.begin();
        std::advance(it, startNodeIndex);
        std::string currentNode = it->first;
        
        // 存储已经访问的边
        std::map<std::string, std::map<std::string, bool>> visitedEdges;
        
        // 存储游走路径
        std::vector<std::string> walkPath;
        walkPath.push_back(currentNode);
        
        // 在第一步游走前就询问用户
        std::string input;
        std::cout << "当前路径: " << currentNode << std::endl;
        std::cout << "按Enter键继续，输入'q'停止游走: ";
        std::getline(std::cin, input);
        if (!input.empty() && (input[0] == 'q' || input[0] == 'Q')) {
            std::cout << "用户停止游走" << std::endl;
            return walkPath;
        }
        while (true) {
            // 检查当前节点是否有出边
            if (adjacencyList[currentNode].empty()) {
                std::cout << "节点 " << currentNode << " 没有出边，游走结束" << std::endl;
                break;
            }
            
            // 随机选择一条出边
            int edgeIndex = std::rand() % adjacencyList[currentNode].size();
            std::string nextNode = adjacencyList[currentNode][edgeIndex].to;
            
            // 检查边是否已经访问过
            if (visitedEdges[currentNode][nextNode]) {
                std::cout << "边 " << currentNode << "->" << nextNode << " 已经访问过，游走结束" << std::endl;
                break;
            }
            
            // 标记边为已访问
            visitedEdges[currentNode][nextNode] = true;
            
            // 更新当前节点，添加到路径
            currentNode = nextNode;
            walkPath.push_back(currentNode);
            
            // 显示当前路径
            std::cout << "当前路径: ";
            for (const auto& node : walkPath) {
                std::cout << node << " ";
            }
            std::cout << std::endl;
            
            // 检查用户输入，决定是否继续
            std::cout << "按Enter键继续，输入'q'停止游走: ";
            std::string input;
            std::getline(std::cin, input);
            if (!input.empty() && (input[0] == 'q' || input[0] == 'Q')) {
                std::cout << "用户停止游走" << std::endl;
                break;
            }
        }
        
        return walkPath;
    }
};



// 在头文件包含后、主函数前添加
bool loadTextFile(const std::string& filename, DirectedGraph& graph);
int main(int argc, char* argv[]) {
    std::string filename;
    DirectedGraph graph;
    bool fileLoaded = false;
    std::string baseName;
    
    // 如果有命令行参数，直接加载文件
    if (argc > 1) {
        filename = argv[1];
        fileLoaded = loadTextFile(filename, graph);
        if (fileLoaded) {
            // 获取文件名（不包括扩展名）作为输出基础名
            baseName = filename;
            size_t lastDot = baseName.find_last_of(".");
            if (lastDot != std::string::npos) {
                baseName = baseName.substr(0, lastDot);
            }
        }
    }
    
    int choice = 0;
    bool running = true;
    
    while (running) {
        // 更新菜单选项
        std::cout << "\n===== 文本图分析工具 =====\n";
        std::cout << "0. 退出程序\n";
        std::cout << "1. 加载文本文件\n";
        std::cout << "2. 显示有向图结构\n";
        std::cout << "3. 导出DOT文件\n";
        std::cout << "4. 生成图像文件\n";
        std::cout << "5. 查询桥接词\n";
        std::cout << "6. 生成新文本\n";
        std::cout << "7. 计算最短路径\n";
        std::cout << "8. 计算PageRank值\n";
        std::cout << "9. 随机游走\n";
        std::cout << "请选择功能 (0-9): ";
        // 处理无效输入
        if (!(std::cin >> choice)) {
            std::cin.clear();  // 清除错误标志
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // 忽略无效输入
            choice = 0;  // 设置无效值
        }
        
        switch (choice) {
            // 保持现有case...
            
            case 1:  // 加载文本文件
                std::cout << "请输入文件名: ";
                std::cin >> filename;
                fileLoaded = loadTextFile(filename, graph);
                if (fileLoaded) {
                    // 获取文件名（不包括扩展名）作为输出基础名
                    baseName = filename;
                    size_t lastDot = baseName.find_last_of(".");
                    if (lastDot != std::string::npos) {
                        baseName = baseName.substr(0, lastDot);
                    }
                    std::cout << "文件加载成功！\n";
                }
                break;
                
            case 2:  // 显示有向图结构
                if (fileLoaded) {
                    graph.printGraph();
                } else {
                    std::cout << "请先加载文本文件！\n";
                }
                break;
                
            case 3:  // 导出DOT文件
                if (fileLoaded) {
                    std::string dotFilename;
                    std::cout << "请输入DOT文件名 [" << baseName << ".dot]: ";
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::getline(std::cin, dotFilename);
                    
                    if (dotFilename.empty()) {
                        dotFilename = baseName + ".dot";
                    }
                    
                    if (graph.exportToDOT(dotFilename)) {
                        std::cout << "DOT文件已生成: " << dotFilename << std::endl;
                    }
                } else {
                    std::cout << "请先加载文本文件！\n";
                }
                break;
                
            case 4:  // 生成图像文件
                if (fileLoaded) {
                    std::string dotFilename = baseName + ".dot";
                    std::string imageFilename;
                    
                    // 先检查DOT文件是否存在
                    std::ifstream checkFile(dotFilename);
                    bool dotExists = checkFile.good();
                    checkFile.close();
                    
                    if (!dotExists) {
                        std::cout << "DOT文件不存在，正在生成...\n";
                        if (!graph.exportToDOT(dotFilename)) {
                            std::cout << "无法生成DOT文件，请先选择'导出DOT文件'功能！\n";
                            break;
                        }
                    }
                    
                    std::cout << "请输入图像文件名 [" << baseName << ".png]: ";
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::getline(std::cin, imageFilename);
                    
                    if (imageFilename.empty()) {
                        imageFilename = baseName + ".png";
                    }
                    
                    if (graph.generateImage(dotFilename, imageFilename)) {
                        std::cout << "可以使用图像查看器打开图形文件: " << imageFilename << std::endl;
                    }
                } else {
                    std::cout << "请先加载文本文件！\n";
                }
                break;
            case 5:  // 查询桥接词
                if (fileLoaded) {
                    std::string word1, word2;
                    std::cout << "请输入第一个单词: ";
                    std::cin >> word1;
                    std::cout << "请输入第二个单词: ";
                    std::cin >> word2;
                    
                    // 转换为小写
                    word1 = toLower(word1);
                    word2 = toLower(word2);
                    
                    // 查找桥接词
                    if (!graph.containsWord(word1) || !graph.containsWord(word2)) {
                        std::cout << "No \"" << word1 << "\" or \"" << word2 << "\" in the graph!" << std::endl;
                    } else {
                        std::vector<std::string> bridges = graph.findBridgeWords(word1, word2);
                        
                        if (bridges.empty()) {
                            std::cout << "No bridge words from \"" << word1 << "\" to \"" << word2 << "\"!" << std::endl;
                        } else {
                            std::cout << "The bridge words from \"" << word1 << "\" to \"" << word2 << "\" are: ";
                            
                            if (bridges.size() == 1) {
                                std::cout << bridges[0] << "." << std::endl;
                            } else {
                                for (size_t i = 0; i < bridges.size() - 1; ++i) {
                                    std::cout << bridges[i];
                                    if (i < bridges.size() - 2) {
                                        std::cout << ", ";
                                    } else {
                                        std::cout << " and ";
                                    }
                                }
                                std::cout << bridges.back() << "." << std::endl;
                            }
                        }
                    }
                } else {
                    std::cout << "请先加载文本文件！\n";
                }
                break;

            case 6:  // 生成新文本
                if (fileLoaded) {
                    std::string inputText;
                    std::cout << "请输入一行新文本: ";
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::getline(std::cin, inputText);
                    
                    std::string newText = graph.generateNewText(inputText);
                    std::cout << "生成的新文本: " << newText << std::endl;
                } else {
                    std::cout << "请先加载文本文件！\n";
                }
                break;
            
            case 7:  // 计算最短路径
                if (fileLoaded) {
                    std::string word1, word2;
                    std::cout << "请输入第一个单词 (源点): ";
                    std::cin >> word1;
                    word1 = toLower(word1);
                    
                    if (!graph.containsWord(word1)) {
                        std::cout << "单词 \"" << word1 << "\" 不在图中！" << std::endl;
                        break;
                    }
                    
                    std::cout << "请输入第二个单词 (终点)，或按Enter键查找到所有单词的最短路径: ";
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::getline(std::cin, word2);
                    
                    if (word2.empty()) {
                        // 计算到所有单词的最短路径
                        std::cout << "计算从 \"" << word1 << "\" 到所有单词的最短路径..." << std::endl;
                        
                        std::map<std::string, PathInfo> allPaths = graph.findAllShortestPaths(word1);
                        
                        if (allPaths.empty()) {
                            std::cout << "没有找到从 \"" << word1 << "\" 出发的路径！" << std::endl;
                        } else {
                            std::cout << "\n从 \"" << word1 << "\" 出发的最短路径:" << std::endl;
                            std::cout << std::string(60, '-') << std::endl;
                            std::cout << std::setw(15) << "目标单词" << " | " 
                                     << std::setw(10) << "路径长度" << " | " 
                                     << "路径" << std::endl;
                            std::cout << std::string(60, '-') << std::endl;
                            
                            for (const auto& path : allPaths) {
                                std::cout << std::setw(15) << path.first << " | " 
                                         << std::setw(10) << path.second.distance << " | ";
                                
                                // 打印路径
                                for (size_t i = 0; i < path.second.path.size(); ++i) {
                                    std::cout << path.second.path[i];
                                    if (i < path.second.path.size() - 1) {
                                        std::cout << "->";
                                    }
                                }
                                std::cout << std::endl;
                            }
                        }
                    } else {
                        // 计算指定两个单词之间的最短路径
                        word2 = toLower(word2);
                        
                        if (!graph.containsWord(word2)) {
                            std::cout << "单词 \"" << word2 << "\" 不在图中！" << std::endl;
                            break;
                        }
                        
                        std::cout << "计算从 \"" << word1 << "\" 到 \"" << word2 << "\" 的最短路径..." << std::endl;
                        
                        PathInfo pathInfo = graph.findShortestPath(word1, word2);
                        
                        if (pathInfo.distance == INT_MAX) {
                            std::cout << "从 \"" << word1 << "\" 到 \"" << word2 << "\" 不可达！" << std::endl;
                        } else {
                            std::cout << "最短路径长度: " << pathInfo.distance << std::endl;
                            std::cout << "路径: ";
                            for (size_t i = 0; i < pathInfo.path.size(); ++i) {
                                std::cout << pathInfo.path[i];
                                if (i < pathInfo.path.size() - 1) {
                                    std::cout << " -> ";
                                }
                            }
                            std::cout << std::endl;
                            
                            // 生成带有高亮路径的图像
                            std::string pathDotFilename = baseName + "_path.dot";
                            std::string pathImageFilename = baseName + "_path.png";
                            
                            if (graph.exportToDOTWithPath(pathDotFilename, pathInfo.path)) {
                                std::cout << "带高亮路径的DOT文件已生成: " << pathDotFilename << std::endl;
                                
                                if (graph.generateImage(pathDotFilename, pathImageFilename)) {
                                    std::cout << "带高亮路径的图像已生成: " << pathImageFilename << std::endl;
                                    std::cout << "可以使用图像查看器打开查看路径。" << std::endl;
                                }
                            }
                        }
                    }
                } else {
                    std::cout << "请先加载文本文件！\n";
                }
                break;
        
            case 8:  // 计算PageRank值
                if (fileLoaded) {
                    std::cout << "PageRank计算选项：\n";
                    std::cout << "1. 标准PageRank\n";
                    std::cout << "2. 改进的PageRank (使用TF-IDF初始化)\n";
                    std::cout << "请选择 (1-2): ";
                    
                    int prChoice;
                    std::cin >> prChoice;
                    
                    std::map<std::string, double> pageRank;
                    
                    // 默认选1，处理空输入
                    if (std::cin.fail()) {
                        std::cin.clear();  // 清除错误标志
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // 忽略无效输入
                        prChoice = 1;  // 默认选择1
                    }

                    
                    if (prChoice == 1) {
                        double dampingFactor;
                        std::cout << "请输入阻尼系数 (推荐0.85): ";
                        std::cin >> dampingFactor;

                        if (dampingFactor <= 0 || dampingFactor >= 1) {
                            std::cout << "阻尼系数必须在0到1之间，将使用默认值0.85\n";
                            dampingFactor = 0.85;
                        }
                        
                        pageRank = graph.calculatePageRank(dampingFactor);
                    } else if (prChoice == 2) {
                        double dampingFactor;
                        std::cout << "请输入阻尼系数 (推荐0.85): ";
                        std::cin >> dampingFactor;

                        if (dampingFactor <= 0 || dampingFactor >= 1) {
                            std::cout << "阻尼系数必须在0到1之间，将使用默认值0.85\n";
                            dampingFactor = 0.85;
                        }
                        
                        pageRank = graph.calculatePageRankWithTFIDF(dampingFactor);
                    } else {
                        std::cout << "无效选择，使用标准PageRank\n";
                        pageRank = graph.calculatePageRank();
                    }
                    
                    // 将结果转换为向量以便排序
                    std::vector<std::pair<std::string, double>> rankedWords;
                    for (const auto& pr : pageRank) {
                        rankedWords.push_back(pr);
                    }
                    
                    // 按PageRank值降序排序
                    std::sort(rankedWords.begin(), rankedWords.end(),
                            [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
                                return a.second > b.second;
                            });
                    
                    // 显示PageRank结果
                    std::cout << "\nPageRank排名结果 (前20):\n";
                    std::cout << std::string(40, '-') << std::endl;
                    std::cout << std::setw(15) << "单词" << " | " 
                            << std::setw(20) << "PageRank值" << std::endl;
                    std::cout << std::string(40, '-') << std::endl;
                    
                    int count = 0;
                    for (const auto& word : rankedWords) {
                        std::cout << std::setw(15) << word.first << " | " 
                                << std::setw(20) << std::fixed << std::setprecision(6) << word.second << std::endl;
                        
                        count++;
                        if (count >= 20) break; // 只显示前20个结果
                    }
                    
                    // 导出PageRank结果到文件
                    std::string prFilename = baseName + "_pagerank.txt";
                    std::ofstream prFile(prFilename);
                    
                    if (prFile.is_open()) {
                        prFile << "单词,PageRank值\n";
                        for (const auto& word : rankedWords) {
                            prFile << word.first << "," << std::fixed << std::setprecision(6) << word.second << "\n";
                        }
                        prFile.close();
                        std::cout << "\nPageRank结果已保存到文件: " << prFilename << std::endl;
                    } else {
                        std::cerr << "无法创建PageRank结果文件\n";
                    }
                    
                } else {
                    std::cout << "请先加载文本文件！\n";
                }
                break;
                // 在switch-case中添加随机游走选项
            case 9:  // 随机游走
                if (fileLoaded) {
                    std::vector<std::string> walkPath = graph.randomWalk();
                    
                    if (!walkPath.empty()) {
                        // 输出游走路径到屏幕
                        std::cout << "\n随机游走完整路径:\n";
                        for (size_t i = 0; i < walkPath.size(); ++i) {
                            std::cout << walkPath[i];
                            if (i < walkPath.size() - 1) {
                                std::cout << " ";
                            }
                        }
                        std::cout << std::endl;
                        
                        // 保存到文件
                        std::string walkFilename = baseName + "_walk.txt";
                        std::ofstream walkFile(walkFilename);
                        
                        if (walkFile.is_open()) {
                            for (size_t i = 0; i < walkPath.size(); ++i) {
                                walkFile << walkPath[i];
                                if (i < walkPath.size() - 1) {
                                    walkFile << " ";
                                }
                            }
                            walkFile.close();
                            std::cout << "随机游走路径已保存到文件: " << walkFilename << std::endl;
                        } else {
                            std::cerr << "无法创建随机游走结果文件\n";
                        }
                    } else {
                        std::cout << "随机游走未产生有效路径\n";
                    }
                } else {
                    std::cout << "请先加载文本文件！\n";
                }
                break;
            case 0:  // 退出程序
                std::cout << "感谢使用！\n";
                running = false;
                break;
                
            // 更新default处理范围
            default:
                std::cout << "无效选择，请重新输入 (0-9)。\n";
                break;
        }
    }
    
    return 0;
}

// 添加加载文本文件的函数
bool loadTextFile(const std::string& filename, DirectedGraph& graph) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return false;
    }
    
    std::string previousWord = "";
    std::string currentWord = "";
    char ch;
    
    while (file.get(ch)) {
        if (std::isalpha(ch)) {
            // 字母直接添加到当前单词
            currentWord += ch;
        } else {
            // 非字母字符处理为空格
            if (!currentWord.empty()) {
                // 单词结束，处理边
                std::string lowerCurrentWord = toLower(currentWord);
                
                if (!previousWord.empty()) {
                    // 添加从前一个单词到当前单词的边
                    graph.addEdge(previousWord, lowerCurrentWord);
                }
                
                previousWord = lowerCurrentWord;
                currentWord = "";
            }
        }
    }
    
    // 处理最后一个单词
    if (!currentWord.empty() && !previousWord.empty()) {
        std::string lowerCurrentWord = toLower(currentWord);
        graph.addEdge(previousWord, lowerCurrentWord);
    }
    
    file.close();
    return true;
}