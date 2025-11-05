#ifndef CHEMS_FUZZY_MAP_HPP
#define CHEMS_FUZZY_MAP_HPP

#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <tuple>
#include <iostream>

#include "trie.hpp"


namespace chm
{
    template <typename value_t>
    class FuzzyMap
    {
    private:
        Trie<value_t, char> trie;
        
        std::vector<std::vector<int>> lev_matrix;
        int max_key_size{0};
        
        const int query_size_multiplier{100};
        const int min_max_cost_thr{8};
        
        int add_cost_middle, add_cost_sides;
        int remove_cost_middle, remove_cost_sides;

    public:

        FuzzyMap(double add_cost_middle = 3.0, double add_cost_sides = 1.0, double remove_cost_middle = 4.0, double remove_cost_sides = 4.0);
        
        void insert_entries(const std::vector<std::pair<std::string, value_t>>& entries);

        std::vector<std::pair<value_t, double>> search(std::string query);
    };


    template<typename value_t>
    FuzzyMap<value_t>::FuzzyMap(double add_cost_middle, double add_cost_sides, double remove_cost_middle, double remove_cost_sides) : trie('\0')
    {
        if(add_cost_middle < 0 || add_cost_sides < 0 || remove_cost_middle < 0 || remove_cost_sides < 0)
            throw std::invalid_argument("Cost must be positive");
        
        this->add_cost_middle = add_cost_middle * this->query_size_multiplier;
        this->add_cost_sides = add_cost_sides * this->query_size_multiplier;

        this->remove_cost_middle = remove_cost_middle * this->query_size_multiplier;
        this->remove_cost_sides = remove_cost_sides * this->query_size_multiplier;

        this->lev_matrix.insert(this->lev_matrix.end(), 2, std::vector<int>(this->max_key_size+1, -1));
    }

    template<typename value_t>
    void FuzzyMap<value_t>::insert_entries(const std::vector<std::pair<std::string, value_t>>& entries)
    {
        int prev_max_key_size = this->max_key_size;
        for (const auto& e : entries)
        {
            this->trie.insert_string(e.first, e.second);
            this->max_key_size = std::max<int>(this->max_key_size, (int)e.first.size());
        }
        
        if(this->max_key_size > prev_max_key_size)
        {
            this->lev_matrix[0].resize(this->max_key_size+1);
            this->lev_matrix[1].resize(this->max_key_size+1);
        }
    }


    static inline std::string to_lower(const std::string& s)
    {
        std::string result = s;
        std::transform(result.begin(), result.end(), result.begin(),
            [](unsigned char c) { return std::tolower(c); });
        
        return result;
    }


    template<typename value_t>
    std::vector<std::pair<value_t, double>> FuzzyMap<value_t>::search(std::string query)
    {
        query = to_lower(query);
        int query_size = query.size();
        
        int max_cost_thr = std::max(this->min_max_cost_thr * this->query_size_multiplier, query_size * this->query_size_multiplier);
        std::unordered_map<value_t, int> value_cost;

        using active_state_t = std::tuple<const Trie<value_t, char>*, int, int>;
        std::stack<active_state_t> active_states;

        struct ActiveStateHash
        {
            std::size_t operator()(const active_state_t& t) const
            {
                auto h1 = std::hash<const Trie<value_t, char>*>{}(std::get<0>(t));
                auto h2 = std::hash<int>{}(std::get<1>(t));
                auto h3 = std::hash<int>{}(std::get<2>(t));
                
                return h1 ^ (h2 << 1) ^ (h3 << 2);
            }
        };

        std::unordered_set<active_state_t, ActiveStateHash> visited_states;

        auto check_add_active_state =
            [&](const Trie<value_t, char> *trie, int pos, int cost)
            {
                auto state = std::make_tuple(trie, pos, cost);
                if(cost <= max_cost_thr && visited_states.find(state) == visited_states.end())
                {
                    active_states.push(state);
                    visited_states.insert(state);
                }
            };

        check_add_active_state(&this->trie, 0, 0.0);
        
        int remove_cost, add_cost;
        while(active_states.size())
        {
            auto [curr_trie, curr_pos, curr_cost] = active_states.top();
            active_states.pop();

            bool is_end_pos = curr_pos == query_size;

            if(curr_pos == 0 || is_end_pos)
            {
                add_cost = this->add_cost_sides;
                remove_cost = this->remove_cost_sides;
            }
            else
            {
                add_cost = this->add_cost_middle;
                remove_cost = this->remove_cost_middle;
            }

            if(is_end_pos)
            {
                for(value_t value : curr_trie->values)
                {
                    auto it = value_cost.find(value);
                    if(it != value_cost.end())
                        value_cost[value] = std::min(value_cost[value], curr_cost);
                    else
                        value_cost[value] = curr_cost;
                }
                
                for(const auto& child_trie : curr_trie->children)
                    check_add_active_state(child_trie, curr_pos, curr_cost + add_cost);
            }
            else
            {
                for(const auto& child_trie : curr_trie->children)
                {
                    if(query[curr_pos] == child_trie->get_key())
                        check_add_active_state(child_trie, curr_pos+1, curr_cost);
                    check_add_active_state(child_trie, curr_pos, curr_cost + add_cost);
                }

                check_add_active_state(curr_trie, curr_pos+1, curr_cost + remove_cost);
            }
        }

        std::vector<std::pair<value_t, double>> result;
        result.reserve(value_cost.size());
        for(const auto& entry : value_cost)
        {
            value_t value = entry.first;
            double cost = (double)entry.second / query_size / this->query_size_multiplier;
            result.push_back(std::make_pair(value, cost));
        }

        std::sort(result.begin(), result.end(), [](const std::pair<value_t, double>& a, const std::pair<value_t, double>& b) { return a.second < b.second; });

        return result;
    }
}

#endif // CHEMS_FUZZY_MAP_HPP