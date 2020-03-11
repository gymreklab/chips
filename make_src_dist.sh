#!/bin/sh

############## Follow these instructions to make a new release ########
# Pull changes
# > git pull origin master
#
# Update the version tag (obviously increment this. use "git tag" to see existing versions)
# > git tag -a v2.0 -m"v2.0"
#
# Push tag update
# > git push origin --tags
#
# Run this to make sure tag gets properly reset
# > ./reconf
#
# Package
# > ./make_src_dist.sh
#
# This will create a tarball file: chips-X.X.XXXX.tar.gz
# Upload this file to Github releases page: https://github.com/gymreklab/chips/releases

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

./configure || die "./configure failed"

GITID=$(./config/git-version-gen .tarball-version) || die "Failed to get GIT ID"
TARBALL=chips-${GITID}.tar.gz
rm -f "$TARBALL"

make distcheck || die "make distcheck failed"
[ -e "$TARBALL" ] || die "can't find Tarball file '$TARBALL' after 'make distcheck'"

echo
echo "Version $GITID is ready for distribution"
echo ""
echo "Upload the following file to GitHub:"
echo "  $TARBALL"
echo ""
echo ""
