#ifndef CHEMS_MISC_HPP
#define CHEMS_MISC_HPP

#include <string>
#include <charconv>
#include <stdexcept>
#include <type_traits>
#include <sstream>


namespace chm
{
    template<typename>
    inline constexpr bool always_false_v = false;


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
        else if (ptr != s.data() + s.size())
        {
            throw std::invalid_argument("Invalid number: " + s);
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


    template<typename T>
    std::string type_to_sql_str(const T& entry)
    {
        if constexpr(std::is_arithmetic_v<T>)
        {
            return numeric_to_str(entry);
        }
        else if constexpr(std::is_same_v<T, std::string>)
        {
            std::string result;
            for(char c : entry)
            {
                if(c == '\\' || c == '"')
                    result += "\\";

                result += c;
            }
            return '"' + result + '"';
        }
        else
            static_assert(always_false_v<T>, "Function's not implemented for a provided type");
    }


    template<typename Iterable, typename CidsBlacklist>
    std::string entries_to_sql_list(const Iterable& entries, const CidsBlacklist& blacklist)
    {
        std::ostringstream cids_stream;
        cids_stream << "{";
        size_t added = 0;
        for(const auto& entry : entries)
        {
            if(blacklist.find(entry) != blacklist.end())
                continue;

            if(added != 0)
                cids_stream << ',';
            cids_stream << type_to_sql_str(entry);
            ++added;
        }
        cids_stream << "}";

        return cids_stream.str();
    }
}

#endif // CHEMS_MISC_HPP
