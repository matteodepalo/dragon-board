# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to, avoiding a bug in XCode 1.5
all.Debug: \
	/Users/matteodepalo/Documents/Università/dragon-board/lib/Debug/libdragon-board.so

all.Release: \
	/Users/matteodepalo/Documents/Università/dragon-board/lib/Release/libdragon-board.so

all.MinSizeRel: \
	/Users/matteodepalo/Documents/Università/dragon-board/lib/MinSizeRel/libdragon-board.so

all.RelWithDebInfo: \
	/Users/matteodepalo/Documents/Università/dragon-board/lib/RelWithDebInfo/libdragon-board.so

# For each target create a dummy rule so the target does not have to exist


# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
/Users/matteodepalo/Documents/Università/dragon-board/lib/Debug/libdragon-board.so:
	/bin/rm -f /Users/matteodepalo/Documents/Università/dragon-board/lib/Debug/libdragon-board.so


/Users/matteodepalo/Documents/Università/dragon-board/lib/Release/libdragon-board.so:
	/bin/rm -f /Users/matteodepalo/Documents/Università/dragon-board/lib/Release/libdragon-board.so


/Users/matteodepalo/Documents/Università/dragon-board/lib/MinSizeRel/libdragon-board.so:
	/bin/rm -f /Users/matteodepalo/Documents/Università/dragon-board/lib/MinSizeRel/libdragon-board.so


/Users/matteodepalo/Documents/Università/dragon-board/lib/RelWithDebInfo/libdragon-board.so:
	/bin/rm -f /Users/matteodepalo/Documents/Università/dragon-board/lib/RelWithDebInfo/libdragon-board.so


