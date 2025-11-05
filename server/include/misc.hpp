#ifndef CHEMS_MISC_HPP
#define CHEMS_MISC_HPP

#include <string>
#include <charconv>
#include <stdexcept>
#include <type_traits>


namespace chm
{
    template<typename T>
    T str_to_numeric(const std::string& s)
    {
        static_assert(std::is_arithmetic_v<T>, "T must be numeric");

        T value{};
        auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);

        if (ec == std::errc::invalid_argument)
        {
            throw std::invalid_argument("Invalid number: " + s);
        }
        else if (ec == std::errc::result_out_of_range)
        {
            throw std::out_of_range("Number out of range: " + s);
        }

        return value;
    }

    template<typename T>
    std::string numeric_to_str(T value) {
        static_assert(std::is_arithmetic_v<T>, "T must be numeric");

        char buffer[64];
        auto [ptr, ec] = std::to_chars(buffer, buffer + sizeof(buffer), value);

        if (ec == std::errc::value_too_large) {
            throw std::runtime_error("Number too large to convert to string");
        }

        return std::string(buffer, ptr);
    }
}

#endif // CHEMS_MISC_HPP