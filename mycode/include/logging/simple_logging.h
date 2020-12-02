#ifndef UTILS_SIMPLE_LOGGING_H
#define UTILS_SIMPLE_LOGGING_H

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>

#define CHECK(x)                                           \
  if (!(x))                                                \
    LogMessageFatal(__FILE__, __LINE__).stream() << "Check failed: " #x << ' '
#define CHECK_LT(x, y) CHECK((x) < (y))
#define CHECK_GT(x, y) CHECK((x) > (y))
#define CHECK_LE(x, y) CHECK((x) <= (y))
#define CHECK_GE(x, y) CHECK((x) >= (y))
#define CHECK_EQ(x, y) CHECK((x) == (y))
#define CHECK_NE(x, y) CHECK((x) != (y))

#ifdef NDEBUG
#define DCHECK(x) ((void)0)
#define DCHECK_LT(x, y) ((void)0)
#define DCHECK_GT(x, y) ((void)0)
#define DCHECK_LE(x, y) ((void)0)
#define DCHECK_GE(x, y) ((void)0)
#define DCHECK_EQ(x, y) ((void)0)
#define DCHECK_NE(x, y) ((void)0)
#else
#define DCHECK(x) CHECK(x)
#define DCHECK_LT(x, y) CHECK((x) < (y))
#define DCHECK_GT(x, y) CHECK((x) > (y))
#define DCHECK_LE(x, y) CHECK((x) <= (y))
#define DCHECK_GE(x, y) CHECK((x) >= (y))
#define DCHECK_EQ(x, y) CHECK((x) == (y))
#define DCHECK_NE(x, y) CHECK((x) != (y))
#endif  // NDEBUG

#define LOG_INFO LogMessage(__FILE__, __LINE__, "Info")
#define LOG_ERROR LogMessage(__FILE__, __LINE__, "Error")
#define LOG_WARNING LogMessage(__FILE__, __LINE__, "Warning")
#define LOG_FATAL LogMessageFatal(__FILE__, __LINE__)
#define LOG(severity) LOG_##severity.stream()

class DateLogger {
public:
    const char *HumanDate() {
        time_t time_value = time(NULL);
        tm now;
        localtime_r(&time_value, &now);
        snprintf(buffer_, sizeof(buffer_), "%02d:%02d:%02d", now.tm_hour, now.tm_min, now.tm_sec);
        return buffer_;
    }

private:
    char buffer_[9];
};

class LogMessage {
public:
    LogMessage(const char *file, int line, const char *level) : log_stream_(std::cerr) {
        log_stream_ << "[" << level << " " << pretty_date_.HumanDate() << "] " << file << ", "
                    << line << ": ";
    }

    virtual ~LogMessage() { log_stream_ << "\n"; }

    std::ostream &stream() { return log_stream_; }

protected:
    std::ostream &log_stream_;

private:
    DateLogger pretty_date_;

    LogMessage(const LogMessage &);

    void operator=(const LogMessage &);
};

class LogMessageFatal : public LogMessage {
public:
    LogMessageFatal(const char *file, int line) : LogMessage(file, line, "Fatal") {}

    ~LogMessageFatal() {
        log_stream_ << "\n";
        abort();
    }

private:
    LogMessageFatal(const LogMessageFatal &);

    void operator=(const LogMessageFatal &);
};

#endif  // SIMPLE_LOGGING_H
