#include "exprtk.hpp"

struct Expression {
    explicit Expression(const std::string &expressionString) {
        symbol_table.add_variable("x",x);
        symbol_table.add_variable("y",y);
        expression.register_symbol_table(symbol_table);
        parser.compile(expressionString, expression);
    }

    double operator()(double _x, double _y) {
        x = _x;
        y = _y;
        auto result = expression.value();
        return result;
    }

private:
    double x{}, y{};
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;
    exprtk::parser<double> parser;
};