#include "../include/vba_interface.h"
#include "../include/pricers.h"
#include "../include/option.h"
#include "../include/pde.h"
#include "../include/fdm.h"
#include "../include/greeks.h"
#include <fstream>

Input get_params_from_file(std::string filename) {
    std::fstream params_stream(filename, std::ios::in);

    double s, t, K, r, v,M,N;
    std::string option_type;

    params_stream >> t >> K >> M >> N >> s >> r >> v >> option_type;    
    params_stream.close();

    Input input;

    input.spot = s
    , input.time_to_maturity = t
    , input.strike = K
    , input.rate = r
    , input.volatility = v
    , input.M = M
    , input.N = N
    , input.option_type = option_type
    ;

    return input;
}

void store_greeks(Output &output,Input input){
    price_inputs in;
    in.spot = input.spot, in.time_to_maturity = input.time_to_maturity, in.strike = input.strike, in.rate = input.rate, in.volatility = input.volatility;
    UniformDiscretization disc = default_UniformDiscretization(in, input.M, input.N);
    in.disc = &disc;

    std::string output_pde = input.option_type + "_greek";
    
    output.delta=delta_option(in,output_pde , input.option_type.c_str());
    output.gamma=gamma_option(in,output_pde , input.option_type.c_str());
    output.theta=theta_option(in,output_pde , input.option_type.c_str());
    output.rho=rho_option(in,output_pde , input.option_type.c_str());
    output.vega=vega_option(in,output_pde , input.option_type.c_str());
}


void create_output_file(std::string filename, Output output) {
    std::remove("output");
    std::remove("output.csv");
    std::ofstream put_out(filename);
    put_out << output.price << '\n';
    put_out << output.delta << '\n';
    put_out << output.gamma << '\n';
    put_out << output.theta << '\n';
    put_out << output.rho << '\n';
    put_out << output.vega << '\n';
    put_out.close();
}

double price_option(Input input) {

    price_inputs in;
    in.spot = input.spot, in.time_to_maturity = input.time_to_maturity, in.strike = input.strike, in.rate = input.rate, in.volatility = input.volatility;
    UniformDiscretization disc = default_UniformDiscretization(in, input.M, input.N);
    in.disc = &disc;

    if (input.option_type == "european_call") {
        return price_european_call(in);
    }
    if (input.option_type == "american_call") {
        return price_american_call(in);
    }
    if (input.option_type == "european_put") {
        return price_european_put(in);
    }
    if (input.option_type == "american_put") {
        return price_american_put(in);
    }
    
    return -100.0;
}


