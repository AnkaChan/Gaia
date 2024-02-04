#pragma once

#include <MeshFrame/Utility/cxxopts.hpp>
#include "../VersionTracker/VersionTracker.h"

namespace GAIA {
	struct CommandParser
	{
		bool gui = false;
		bool showGitInfo = false;
		bool runOnCPU = false;
		std::string recoveryStateFile = "";
		std::string repoRoot = "";

		CommandParser() :
			options("EBDAppParams", "GAIA physics application.")
		{
			options.add_options()
				("g,gui", "", cxxopts::value<bool>(gui))
				("r,recoveryState", "", cxxopts::value<std::string>(recoveryStateFile))
				("R,repoRoot", "", cxxopts::value<std::string>(repoRoot))
				("CPU", "", cxxopts::value<bool>(runOnCPU))
				("git", "", cxxopts::value<bool>(showGitInfo))
				;
		}

		void parse(int argc, char** argv) {
			try
			{
				auto result = options.parse(argc, argv);
			}
			catch (const cxxopts::OptionException& e)
			{
				std::cout << "error parsing options: " << e.what() << std::endl;
				std::cout << options.help();
				exit(1);
			}
		}

		cxxopts::Options options;
	};
}

