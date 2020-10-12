#include <iostream>
#include <vector>
using namespace std;

vector<int> twoSum(vector<int>& nums, int target)
{
    vector<int> ans;
    int tmp = 0;
    int max_value = 0;
    int length = nums.size();
    for(int i = 0; i < length; i++)
    {
        if(nums[i] > 0)
            tmp = nums[i];
        else
            tmp = -nums[i];
            
        if(tmp >= max_value)
            max_value = tmp;
    }

    vector<int> tmp_num;
    int dim = max_value + 1;
    tmp_num.resize(4 * dim);

    for(int i = 0; i < length; i++)
    {
        int p = nums[i];
        int a = target - p;
        if(a >= 0)
        {
            if (tmp_num[0 * dim + a] == 1)
            {
                ans.push_back(i);
                ans.push_back(tmp_num[2 * dim + a]);
                return ans;
            }
            tmp_num[0 * dim + p] = 1;
            tmp_num[2 * dim + p] = i;

        }
        else
        {
            a = -a;
            if(tmp_num[1 * dim + a] == 1)
            {
                ans.push_back(i);
                ans.push_back(tmp_num[3 * dim + a]);
                return ans;
            }
            tmp_num[1 * dim + p] = 1;
            tmp_num[3 * dim + p] = i;
        }
    }
}

int main()
{
    vector<int> nums = {2, 7, 11, 15};

    int target = 9;

    vector<int> ans = twoSum(nums, target);

    cout << ans[0] << endl;
    cout << ans[1] << endl;

}