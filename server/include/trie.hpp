#ifndef CHEMS_TRIE_HPP
#define CHEMS_TRIE_HPP

#include <memory_resource>
#include <vector>
#include <utility>
#include <string>
#include <set>


namespace chm
{
    template<typename value_t>
    class Trie
    {
    private:
        std::pmr::memory_resource *pool;
        bool own_pool;

        std::string s{""};
        char c;

        typename std::pmr::vector<Trie<value_t>>::iterator insert_char(char c);

    public:

        std::pmr::vector<Trie> children;
        std::set<value_t> values;

        Trie(char c, std::pmr::memory_resource *pool = nullptr);

        ~Trie();

        void insert(std::string key, value_t value);

        std::string get_string() const;
        char get_char() const;

        bool operator<(char b);
    };


    template<typename value_t>
    Trie<value_t>::Trie(char c, std::pmr::memory_resource *pool) : c(c)
    {
        this->own_pool = !pool;
        if (this->own_pool)
            this->pool = new std::pmr::monotonic_buffer_resource();
        else 
            this->pool = pool;

        this->children = std::pmr::vector<Trie>(this->pool);
        this->children.reserve(26);
    }


    template<typename value_t>
    Trie<value_t>::~Trie()
    {
        if(this->own_pool)
            delete this->pool;
    }


    template<typename value_t>
    typename std::pmr::vector<Trie<value_t>>::iterator Trie<value_t>::insert_char(char c)
    {
        auto it = std::lower_bound(this->children.begin(), this->children.end(), c);
        
        if(it == this->children.end() || it->c != c)
        {
            Trie trie = Trie(c, this->pool);
            trie.s = this->s + c;
            it = this->children.insert(it, trie);
        }
        
        return it;
    }


    template<typename value_t>
    void Trie<value_t>::insert(std::string key, value_t value)
    {
        Trie<value_t> *curr_node = this;
        for(char c : key)
            curr_node = &(*curr_node->insert_char(std::tolower(c)));
        curr_node->values.insert(value);
    }


    template<typename value_t>
    std::string Trie<value_t>::get_string() const { return this->s; }

    template<typename value_t>
    char Trie<value_t>::get_char() const { return this->c; }


    template<typename value_t>
    bool Trie<value_t>::operator<(char c) { return this->c < c; }
}

#endif // CHEMS_TRIE_HPP