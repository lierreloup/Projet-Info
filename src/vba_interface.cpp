#include "../include/vba_interface.h"
#include <fstream>

Input get_params_from_file(std::string filename) {
    std::fstream params_stream(filename, std::ios::in);

    double s, t, K, r, v;
    params_stream >> s >> t >> K >> r >> v;    
    params_stream.close();

    Input input;

    input.spot = s
    , input.time_to_maturity = t
    , input.strike = K
    , input.rate = r
    , input.volatility = v
    ;

    return input;
}

void create_output_file(std::string filename, Output output) {
    std::ofstream put_out(filename);
    put_out << output.price << '\n';
    put_out << output.delta << '\n';

    put_out.close();
}