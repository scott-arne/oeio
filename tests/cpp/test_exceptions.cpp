#include <gtest/gtest.h>

#include "oeio/oeio.h"

#include <string>
#include <filesystem>

namespace oeio {
namespace test {

class ExceptionsTest : public ::testing::Test {
protected:
    void SetUp() override {
        tmpdir_ = std::filesystem::temp_directory_path() /
                  ("oeio_exc_test_" + std::to_string(getpid()));
        std::filesystem::create_directories(tmpdir_);
    }

    void TearDown() override {
        std::filesystem::remove_all(tmpdir_);
    }

    std::filesystem::path tmpdir_;
};

TEST_F(ExceptionsTest, UnknownWriteExtensionThrowsFormatError) {
    std::string path = (tmpdir_ / "bad.xyzzy").string();
    EXPECT_THROW(oeio::write(path), oeio::FormatError);
}

TEST_F(ExceptionsTest, UnknownReadExtensionThrowsFormatError) {
    std::string path = (tmpdir_ / "bad.xyzzy").string();
    EXPECT_THROW(oeio::read(path), oeio::FormatError);
}

TEST_F(ExceptionsTest, MissingFileThrowsFileError) {
    std::string path = (tmpdir_ / "does_not_exist.sdf").string();
    EXPECT_THROW(oeio::read(path), oeio::FileError);
}

TEST_F(ExceptionsTest, FormatErrorIsAnError) {
    std::string path = (tmpdir_ / "bad.xyzzy").string();
    try {
        oeio::write(path);
        FAIL() << "Expected FormatError";
    } catch (const oeio::Error&) {
        SUCCEED();
    }
}

TEST_F(ExceptionsTest, ErrorIsAStdException) {
    std::string path = (tmpdir_ / "bad.xyzzy").string();
    try {
        oeio::write(path);
        FAIL() << "Expected exception";
    } catch (const std::exception& e) {
        EXPECT_NE(std::string(e.what()).find("unrecognized"), std::string::npos);
    }
}

}  // namespace test
}  // namespace oeio
