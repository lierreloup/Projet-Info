// Your First C++ Program
#include <fstream>
int main() {
    // Create an output filestream object
    std::ofstream myFile("C:\\Users\\UTILISATEUR\\projects\\Projet-Info\\vba\\foo.csv");
    
    // Send data to the stream
    myFile << "Foo\n";
    myFile << "1\n";
    myFile << "2\n";
    myFile << "3\n";
    
    // Close the file
    myFile.close();
    
    return 0;
}

