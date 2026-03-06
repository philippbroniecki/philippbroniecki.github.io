// @ts-check
import { defineConfig } from 'astro/config';
import sitemap from '@astrojs/sitemap';
import mdx from '@astrojs/mdx';

export default defineConfig({
  site: 'https://philippbroniecki.com',
  integrations: [sitemap(), mdx()],
  redirects: {
    '/research/': '/academia/research/',
    '/teaching/': '/academia/teaching/',
  },
});
