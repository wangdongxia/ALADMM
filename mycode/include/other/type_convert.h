#ifndef UTILS_TYPE_CONVERT_H
#define UTILS_TYPE_CONVERT_H

#include <sstream>
#include <string>

#include "logging/simple_logging.h"

template <typename Target, typename Source, bool Same>
class Converter{
public:
    static Target Convert(const Source &arg){
        Target ret;
        std::stringstream ss;
        if (!(ss<<arg && ss>>ret && ss.eof())) {
            LOG(FATAL) << "类型转换失败";
        }
        return ret;
    }
};

template <typename Target, typename Source>
class Converter<Target, Source, true>{
public:
    static Target Convert(const Source &arg){
        return arg;
    }
};

template <typename Source>
class Converter<std::string, Source, false>{
public:
    static std::string Convert(const Source &arg){
        std::ostringstream ss;
        ss<<arg;
        return ss.str();
    }
};

template <typename Target>
class Converter<Target, std::string, false>{
public:
    static Target Convert(const std::string &arg){
        Target ret;
        std::istringstream ss(arg);
        if (!(ss>>ret && ss.eof())) {
            LOG(FATAL) << "类型转换失败";
        }
        return ret;
    }
};

template <typename T1, typename T2>
struct IsSame {
    static const bool value = false;
};

template <typename T>
struct IsSame<T, T>{
    static const bool value = true;
};

template<typename Target, typename Source>
Target Convert(const Source &arg)
{
    return Converter<Target, Source, IsSame<Target, Source>::value>::Convert(arg);
}

#endif
