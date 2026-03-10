#!/usr/bin/env python3
"""
从LHCb ThOr Functors参考页面提取函数名和描述
输出为CSV文件
"""

import json
import re
import csv
from html.parser import HTMLParser
from html import unescape


class FunctorParser(HTMLParser):
    """解析ThOr Functors页面的HTML解析器"""
    
    def __init__(self):
        super().__init__()
        self.functors = []
        self.current_functor = None
        self.in_sig_name = False
        self.in_description = False
        self.in_dd = False
        self.description_text = ""
        self.functor_name = ""
        self.depth = 0
        self.dd_depth = 0
        
    def handle_starttag(self, tag, attrs):
        attrs_dict = dict(attrs)
        
        # 检测函数名 - 在 dt class="sig sig-object py" 下的 span.sig-name descname
        if tag == 'span':
            classes = attrs_dict.get('class', '').split()
            if 'sig-name' in classes and 'descname' in classes:
                self.in_sig_name = True
                self.functor_name = ""
        
        # 检测描述 - 在 dd 标签内的第一个 p 标签
        if tag == 'dd':
            self.in_dd = True
            self.dd_depth = 0
            self.description_text = ""
        
        if tag == 'p' and self.in_dd and self.dd_depth == 0:
            self.in_description = True
            self.description_text = ""
            
        # 跟踪嵌套深度
        if tag in ('dd', 'dl', 'div'):
            if self.in_dd:
                self.dd_depth += 1
    
    def handle_endtag(self, tag):
        if tag == 'span' and self.in_sig_name:
            self.in_sig_name = False
            # 保存函数名，等待描述
            self.current_functor = self.functor_name.strip()
            
        if tag == 'p' and self.in_description:
            self.in_description = False
            # 描述结束，保存结果
            if self.current_functor and self.description_text.strip():
                # 清理描述文本
                desc = self.clean_description(self.description_text)
                if desc:
                    self.functors.append({
                        'name': self.current_functor,
                        'description': desc
                    })
            self.current_functor = None
            
        if tag == 'dd':
            self.in_dd = False
            self.dd_depth = 0
            
        if tag in ('dd', 'dl', 'div') and self.in_dd:
            self.dd_depth -= 1
    
    def handle_data(self, data):
        if self.in_sig_name:
            self.functor_name += data
        if self.in_description:
            self.description_text += data
    
    def clean_description(self, text):
        """清理描述文本"""
        # 解码HTML实体
        text = unescape(text)
        # 移除多余空白
        text = ' '.join(text.split())
        # 移除特殊字符
        text = text.replace('\xa0', ' ')
        return text.strip()


def parse_html_file(html_content):
    """解析HTML内容并提取函数信息"""
    parser = FunctorParser()
    parser.feed(html_content)
    return parser.functors


def parse_html_with_regex(html_content):
    """使用正则表达式解析HTML（备用方法）"""
    functors = []
    
    # 查找所有函数定义块
    # 模式：<dt ... id="Functors.XXX">...<span class="pre">XXX</span>...</dt>
    # 然后找紧跟的<dd>...<p>描述</p>...
    
    # 匹配函数名
    name_pattern = r'<dt[^>]*id="Functors\.(\w+)"[^>]*>.*?<span class="sig-name descname">.*?<span class="pre">(\w+)</span>.*?</span>'
    
    # 更精确的模式：找到整个函数定义块
    block_pattern = r'<dl class="py function">.*?<dt[^>]*id="Functors\.(\w+)"[^>]*>.*?</dt>\s*<dd>\s*<p>([^<]*(?:<(?!/p>)[^<]*)*)</p>'
    
    # 简化模式：直接匹配函数名和第一个描述段落
    # 函数名模式
    func_name_pattern = r'<span class="sig-name descname">\s*<span class="pre">(\w+)</span>\s*</span>'
    # 找到所有函数名
    func_names = re.findall(func_name_pattern, html_content)
    
    # 描述模式 - 在<dd>后的第一个<p>标签内容
    # 使用更宽松的匹配
    dd_p_pattern = r'<dd>\s*<p>([^<]+(?:<(?!/p>)[^<]*)*)</p>'
    descriptions = re.findall(dd_p_pattern, html_content)
    
    # 合并结果
    for i, name in enumerate(func_names):
        if i < len(descriptions):
            desc = unescape(descriptions[i])
            desc = ' '.join(desc.split())
            desc = desc.replace('\xa0', ' ')
            functors.append({
                'name': name,
                'description': desc
            })
    
    return functors


def parse_html_advanced(html_content):
    """高级HTML解析方法"""
    functors = []
    
    # 分割成函数块
    # 每个函数以 <dl class="py function"> 开始
    blocks = re.split(r'<dl class="py function">', html_content)
    
    for block in blocks[1:]:  # 跳过第一个空块
        # 提取函数名
        name_match = re.search(
            r'<span class="sig-name descname">\s*<span class="pre">(\w+)</span>',
            block
        )
        
        if not name_match:
            continue
            
        name = name_match.group(1)
        
        # 提取描述 - 在<dd>后的第一个<p>标签
        # 找到<dd>标签
        dd_match = re.search(r'<dd>(.*?)(?:</dd>|<dl)', block, re.DOTALL)
        if dd_match:
            dd_content = dd_match.group(1)
            # 找第一个<p>标签
            p_match = re.search(r'<p>([^<]*(?:<(?!/p>)[^<]*)*)</p>', dd_content, re.DOTALL)
            if p_match:
                desc = p_match.group(1)
                # 清理描述
                desc = unescape(desc)
                desc = re.sub(r'<[^>]+>', '', desc)  # 移除HTML标签
                desc = ' '.join(desc.split())  # 清理空白
                desc = desc.replace('\xa0', ' ')
                
                if desc and not desc.startswith('C++ Representation'):
                    functors.append({
                        'name': name,
                        'description': desc
                    })
    
    return functors


def remove_duplicates(functors):
    """移除重复项，保留第一次出现的"""
    seen = set()
    unique = []
    for f in functors:
        if f['name'] not in seen:
            seen.add(f['name'])
            unique.append(f)
    return unique


def main():
    # 读取JSON文件
    json_file = './thor_functors_page.json'
    output_csv = './data/thor_functors.csv'
    
    print(f"读取JSON文件: {json_file}")
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    html_content = data.get('data', {}).get('html', '')
    
    if not html_content:
        print("错误：无法从JSON文件中获取HTML内容")
        return
    
    print("解析HTML内容...")
    
    # 使用多种方法解析
    functors1 = parse_html_advanced(html_content)
    functors2 = parse_html_with_regex(html_content)
    
    # 合并结果，使用第一个非空结果
    if functors1:
        functors = functors1
        print(f"方法1找到 {len(functors1)} 个函数")
    elif functors2:
        functors = functors2
        print(f"方法2找到 {len(functors2)} 个函数")
    else:
        print("警告：两种方法都未找到函数")
        functors = []
    
    # 移除重复
    functors = remove_duplicates(functors)
    print(f"去重后剩余 {len(functors)} 个函数")
    
    # 写入CSV
    print(f"写入CSV文件: {output_csv}")
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=',')  # 使用制表符分隔
        writer.writerow(['Function', 'Discription'])
        for func in functors:
            writer.writerow([func['name'], func['description']])
    
    print(f"完成！共提取 {len(functors)} 个函数")
    
    # 显示前10个作为示例
    print("\n前10个函数示例：")
    print("-" * 80)
    for func in functors[:10]:
        print(f"{func['name']}\t{func['description'][:60]}...")


if __name__ == '__main__':
    main()
