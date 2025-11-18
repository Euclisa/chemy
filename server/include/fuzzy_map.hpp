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
#include <roaring/roaring64map.hh>

#include "trie.hpp"


namespace chm
{
    template <typename value_t>
    class FuzzyMap
    {
    private:
        double similarity_thr;
        size_t qgram_len; // q-gram length

        // Indexed by entry index (0..n-1) -> (bitmap, popcount)
        using bitmap_entry_t = std::pair<roaring::Roaring64Map, size_t>;
        std::vector<bitmap_entry_t> entry_bitmaps;

        // Entry storage: index -> (original_key, value)
        std::vector<std::pair<std::string, value_t>> entries;

        // Reverse lookup: qgram string -> bit position (dense 0,1,2...)
        std::unordered_map<std::string, size_t> qgram_to_bit_idx;


        static inline std::string preprocess(const std::string& s);


        std::pair<roaring::Roaring64Map, size_t> qgram_bitmap_from_str_populate(const std::string& str);
        std::pair<roaring::Roaring64Map, size_t> qgram_bitmap_from_str(const std::string& str);

    public:

        FuzzyMap(const std::vector<std::pair<std::string, value_t>>& entries, double similarity_thr = 0.3, size_t q = 3);

        std::vector<std::pair<value_t, double>> search(std::string query);
    };


    template<typename value_t>
    FuzzyMap<value_t>::FuzzyMap(const std::vector<std::pair<std::string, value_t>>& entries, double similarity_thr, size_t qgram_len) :
    similarity_thr(similarity_thr), qgram_len(qgram_len), entries(entries)
    {
        if(this->qgram_len == 0) throw std::invalid_argument("q must be > 0");

        this->entry_bitmaps.reserve(entries.size());
        std::transform(this->entries.begin(), this->entries.end(), std::back_inserter(this->entry_bitmaps),
            [this](const std::pair<std::string, value_t>& entry)
            { return this->qgram_bitmap_from_str_populate(entry.first); });
    }


    template <typename value_t>
    inline std::string FuzzyMap<value_t>::preprocess(const std::string& s)
    {
        std::string result = s;
        std::transform(result.begin(), result.end(), result.begin(),
            [](unsigned char c) { return std::tolower(c); });
        
        result.erase(std::remove_if(result.begin(), result.end(),
            [](unsigned char c) { return std::isspace(c); }),
            result.end());
        
        return result;
    }


    template<typename value_t>
    std::pair<roaring::Roaring64Map, size_t> FuzzyMap<value_t>::qgram_bitmap_from_str_populate(const std::string& str)
    {
        if(str.size() == 0)
            throw std::invalid_argument("'str' argument must be non-empty string");

        std::string str_processed = this->preprocess(str);
        size_t str_size = str_processed.size();
        size_t upper_bound = str_size < this->qgram_len ? 0 : str_size - this->qgram_len;
        roaring::Roaring64Map bitmap;
        for(size_t i = 0; i <= upper_bound; ++i)
        {
            std::string qgram = str_processed.substr(i, this->qgram_len);
            auto qgram_bit_idx_it = this->qgram_to_bit_idx.find(qgram);
            size_t qgram_bit_idx;
            if(qgram_bit_idx_it == this->qgram_to_bit_idx.end())
            {
                qgram_bit_idx = this->qgram_to_bit_idx.size();
                this->qgram_to_bit_idx[qgram] = qgram_bit_idx;
            }
            else
                qgram_bit_idx = qgram_bit_idx_it->second;

            bitmap.add(qgram_bit_idx);
        }

        return std::make_pair(bitmap, bitmap.cardinality());
    }


    template<typename value_t>
    std::pair<roaring::Roaring64Map, size_t> FuzzyMap<value_t>::qgram_bitmap_from_str(const std::string& str)
    {
        if(str.size() == 0)
            throw std::invalid_argument("'str' argument must be non-empty string");

        std::string str_processed = this->preprocess(str);
        size_t str_size = str_processed.size();
        size_t upper_bound = str_size < this->qgram_len ? 0 : str_size - this->qgram_len + 1;
        roaring::Roaring64Map bitmap;
        for(size_t i = 0; i < upper_bound; ++i)
        {
            std::string qgram = str_processed.substr(i, this->qgram_len);
            auto qgram_bit_idx_it = this->qgram_to_bit_idx.find(qgram);
            if(qgram_bit_idx_it != this->qgram_to_bit_idx.end())
                bitmap.add(qgram_bit_idx_it->second);
        }

        return std::make_pair(bitmap, bitmap.cardinality());
    }


    template<typename value_t>
    std::vector<std::pair<value_t, double>> FuzzyMap<value_t>::search(std::string query)
    {
        auto [query_bitmap, query_popcount] = this->qgram_bitmap_from_str(query);
        std::vector<std::pair<value_t, double>> results;
        std::unordered_set<value_t> accepted_values;
        for(size_t i = 0; i < this->entry_bitmaps.size(); ++i)
        {
            const auto& [key_bitmap, key_popcount] = this->entry_bitmaps[i];
            auto and_bitmap = query_bitmap & key_bitmap;
            size_t and_popcount = and_bitmap.cardinality();
            size_t or_popcount = query_popcount + key_popcount - and_popcount;

            value_t value = this->entries[i].second;
            double similarity = (double)and_popcount / (double)or_popcount; // 'or_popcount' is never zero
            if(similarity >= this->similarity_thr && accepted_values.find(value) == accepted_values.end())
            {
                results.push_back(std::make_pair(this->entries[i].second, similarity));
                accepted_values.insert(value);
            }
        }

        std::sort(results.begin(), results.end(),
            [](const std::pair<value_t, double>& a, const std::pair<value_t, double>& b)
            { return a.second > b.second; });
        
        return results;
    }

}

#endif // CHEMS_FUZZY_MAP_HPP