#ifndef CHEMS_TRIE_HPP
#define CHEMS_TRIE_HPP

#include <memory_resource>
#include <vector>
#include <utility>
#include <string>
#include <set>
#include <algorithm>


namespace chm
{
    template<typename value_t, typename key_t>
    class Trie
    {
    private:

        key_t key;

        Trie *parent;

        bool own_pool;
        std::pmr::memory_resource *pool;
        std::pmr::polymorphic_allocator<Trie> alloc;

        uint32_t depth;

    public:

        std::vector<Trie*> children;
        std::set<value_t> values;

        Trie(key_t key, Trie *parent = nullptr, std::pmr::memory_resource *pool = nullptr);

        ~Trie();

        Trie *insert_key(key_t key);
        
        template<typename T = key_t>
        requires std::is_same_v<T, char>
        void insert_string(std::string key, value_t value);

        template<typename T = key_t>
        requires std::is_same_v<T, char>
        std::string get_string() const;

        std::vector<key_t> get_path(bool include_root = false) const;

        key_t get_key() const;
        uint32_t get_depth() const;
        Trie *get_parent() const;
    };


    template<typename value_t, typename key_t>
    Trie<value_t, key_t>::Trie(key_t key, Trie *parent, std::pmr::memory_resource *pool)
        : key(key),
        parent(parent),
        own_pool(!pool),
        pool(pool ? pool : new std::pmr::monotonic_buffer_resource()),
        alloc(this->pool),
        depth(parent ? parent->depth + 1 : 0) {}


    template<typename value_t, typename key_t>
    Trie<value_t, key_t>::~Trie()
    {
        if(this->own_pool)
            delete this->pool;
    }


    template<typename value_t, typename key_t>
    Trie<value_t, key_t>* Trie<value_t, key_t>::insert_key(key_t key)
    {
        auto it = std::lower_bound(this->children.begin(), this->children.end(), key, [](const Trie<value_t, key_t> *a, key_t b) { return a->get_key() < b; });
        
        if(it == this->children.end() || (*it)->get_key() != key)
        {
            Trie* trie = alloc.template new_object<Trie>(key, this, this->pool);
            it = this->children.insert(it, trie);
        }
        
        return *it;
    }


    template<typename value_t, typename key_t>
    template<typename T>
    requires std::is_same_v<T, char>
    void Trie<value_t, key_t>::insert_string(std::string key, value_t value)
    {
        Trie<value_t, char> *curr_node = this;
        for(char c : key)
            curr_node = curr_node->insert_key(std::tolower(c));
        curr_node->values.insert(value);
    }


    template<typename value_t, typename key_t>
    template<typename T>
    requires std::is_same_v<T, char>
    std::string Trie<value_t, key_t>::get_string() const
    {
        std::vector<char> path = this->get_path();
        return std::string(path.begin(), path.end());
    }


    template<typename value_t, typename key_t>
    std::vector<key_t> Trie<value_t, key_t>::get_path(bool include_root) const
    {
        uint32_t path_len = include_root ? this->depth+1 : this->depth;
        if(path_len == 0)
            return std::vector<key_t>();

        std::vector<key_t> path(path_len, this->key);
        const Trie *curr_node = this;
        uint32_t curr_i = path_len-1;
        while((curr_node && include_root) || (curr_node->parent && !include_root))
        {
            path[curr_i--] = curr_node->key;
            curr_node = curr_node->parent;
        }

        return path;
    }


    template<typename value_t, typename key_t>
    key_t Trie<value_t, key_t>::get_key() const { return this->key; }


    template<typename value_t, typename key_t>
    uint32_t Trie<value_t, key_t>::get_depth() const { return this->depth; }

    template<typename value_t, typename key_t>
    Trie<value_t, key_t> *Trie<value_t, key_t>::get_parent() const { return this->parent; }
}

#endif // CHEMS_TRIE_HPP
