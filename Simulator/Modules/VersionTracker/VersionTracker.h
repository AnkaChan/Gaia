#pragma once

#include "../../3rdParty/cmake-git-version-tracking/git.h"

#include <string>
#include <sstream>
#include <iostream>

struct GitVersionTracker
{
	static std::string getGitInfo() {
		std::stringstream ss;
		//ss << "----------------------------------------------------\n";
		ss << "Git info at the time of compilation:\n"
			;

		ss << " - Git hash: " << git_CommitSHA1() << "\n";
		ss << " - Commit date: " << git_CommitDate() << "\n";
		ss << " - Has uncommitted changes at time of build: " << 
			(git_AnyUncommittedChanges() ? "true" : "false") << "\n";
		ss << " - Commit subject: " << git_CommitSubject() << "\n";
		ss << " - Commit body: " << git_CommitBody() << "\n";
		//std::cout << "Commit description:" << git_Describe() << "\n";
		//ss << "----------------------------------------------------\n";

		return ss.str();
	}

	static void printGitInfo() {
		std::cout << "----------------------------------------------------\n";
		std::cout << getGitInfo();
		std::cout << "----------------------------------------------------\n";

	}
};

